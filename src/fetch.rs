// SPDX-License-Identifier: MIT
//! Reading inputs: local files or URLs (through the system's curl), gzipped or not (through
//! gzip). No crates; both tools are on every macOS and Linux machine.
use std::io::Write;
use std::process::{Command, Stdio};

pub const BKG: &str = "https://igs.bkg.bund.de/root_ftp/IGS";

pub fn fetch(url: &str) -> Result<Vec<u8>, String> {
    let out = Command::new("curl").args(["-sSfL", "--max-time", "180", url]).output()
        .map_err(|e| format!("curl: {e}"))?;
    if !out.status.success() {
        return Err(format!("{url}: {}", String::from_utf8_lossy(&out.stderr).trim()));
    }
    Ok(out.stdout)
}

fn gunzip(data: &[u8]) -> Result<Vec<u8>, String> {
    let mut child = Command::new("gzip").arg("-dc").stdin(Stdio::piped()).stdout(Stdio::piped())
        .spawn().map_err(|e| format!("gzip: {e}"))?;
    let mut stdin = child.stdin.take().unwrap();
    let input = data.to_vec();
    let writer = std::thread::spawn(move || stdin.write_all(&input));
    let out = child.wait_with_output().map_err(|e| format!("gzip: {e}"))?;
    writer.join().unwrap().map_err(|e| format!("gzip: {e}"))?;
    if !out.status.success() { return Err("gzip failed".into()); }
    Ok(out.stdout)
}

/// Text of a local file or URL, decompressed if gzipped.
pub fn read_text(src: &str) -> Result<String, String> {
    let data = if src.contains("://") { fetch(src)? } else {
        std::fs::read(src).map_err(|e| format!("{src}: {e}"))?
    };
    let data = if data.starts_with(&[0x1f, 0x8b]) { gunzip(&data)? } else { data };
    Ok(String::from_utf8_lossy(&data).into_owned())
}

/// URLs of the newest IGS ultra-rapid orbit and its ERP file on BKG.
pub fn latest_igu(now_gps: f64) -> Result<(String, String), String> {
    let week = (now_gps / crate::time::WEEK).floor() as i64;
    for w in [week, week - 1] {
        let listing = match fetch(&format!("{BKG}/products/{w}/")) { Ok(b) => String::from_utf8_lossy(&b).into_owned(), Err(_) => continue };
        let mut names: Vec<&str> = listing.match_indices("IGS0OPSULT_")
            .map(|(i, _)| &listing[i..(i + 41).min(listing.len())])
            .filter(|s| s.len() == 41 && s.ends_with("_02D_15M_ORB.SP3.gz"))
            .collect();
        names.sort();
        names.dedup();
        if let Some(orb) = names.last() {
            let erp = orb.replace("15M_ORB.SP3", "01D_ERP.ERP");
            return Ok((format!("{BKG}/products/{w}/{orb}"), format!("{BKG}/products/{w}/{erp}")));
        }
    }
    Err("no IGS ultra-rapid orbit found on BKG".into())
}

/// Candidate URLs for one day's merged broadcast ephemeris on BKG.
pub fn brdc_urls(y: i64, doy: i64) -> Vec<String> {
    let base = format!("{BKG}/BRDC/{y}/{doy:03}/");
    vec![format!("{base}BRDC00WRD_S_{y}{doy:03}0000_01D_MN.rnx.gz"),
         format!("{base}BRDC00WRD_R_{y}{doy:03}0000_01D_MN.rnx.gz"),
         format!("{base}BRDC00IGS_R_{y}{doy:03}0000_01D_MN.rnx.gz")]
}

pub const CODE: &str = "https://zhw-b.s3.cloud.switch.ch/aiub/CODE";

/// URLs of CODE's newest 0-5 day orbit prediction (GPS, GLONASS, Galileo, every 5 minutes) and its
/// ERP file, from the University of Bern's public bucket (the S3 listing API; the bucket has no
/// index pages, and the old ftp.aiub.unibe.ch no longer answers).
pub fn latest_code_prd(now: f64) -> Result<(String, String), String> {
    let listing = fetch(&format!("https://zhw-b.s3.cloud.switch.ch/aiub?list-type=2&prefix=CODE/COD0OPSPRD_2&max-keys=1000"))?;
    let listing = String::from_utf8_lossy(&listing);
    let mut names: Vec<&str> = listing.split("<Key>CODE/").skip(1)
        .filter_map(|s| s.split('<').next())
        .filter(|n| n.ends_with("_05D_05M_ORB.SP3"))
        .collect();
    names.sort();
    let orb = names.last().ok_or("no CODE 5-day prediction found")?;
    let erp = orb.replace("_05D_05M_ORB.SP3", "_21D_06H_ERP.ERP");
    // The file named for day N starts at 00:00 of day N+1 and is issued during day N, so for most
    // of a day the newest one starts in the future. The one before it then covers now; both go in,
    // the newer taking over from its first epoch.
    let starts = |n: &str| -> Option<f64> {
        let (y, d): (i64, i64) = (n.get(11..15)?.parse().ok()?, n.get(15..18)?.parse().ok()?);
        Some(((crate::time::days_from_civil(y, 1, 1) + d) as f64) * 86400.0 - crate::time::UNIX_GPS0)
    };
    let mut orbits = vec![format!("{CODE}/{orb}")];
    if starts(orb).map_or(false, |t| t > now - 6.0 * 3600.0) && names.len() >= 2 {
        orbits.insert(0, format!("{CODE}/{}", names[names.len() - 2]));
    }
    Ok((orbits.join(","), format!("{CODE}/{erp}")))
}
