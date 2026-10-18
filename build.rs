// SPDX-License-Identifier: MIT
//! Compile the vendored ERFA C sources (vendor/erfa, BSD-3-Clause) into a static library with the
//! system C compiler. No crates: CC (default "cc") and AR (default "ar") must be on the path.
use std::{env, fs, path::PathBuf, process::Command};

fn main() {
    let out = PathBuf::from(env::var("OUT_DIR").unwrap());
    let cc = env::var("CC").unwrap_or_else(|_| "cc".into());
    let ar = env::var("AR").unwrap_or_else(|_| "ar".into());
    let mut objects = Vec::new();
    let mut sources: Vec<PathBuf> = fs::read_dir("vendor/erfa").unwrap()
        .map(|e| e.unwrap().path())
        .filter(|p| p.extension().map_or(false, |x| x == "c"))
        .collect();
    sources.sort();
    for src in &sources {
        let obj = out.join(src.file_stem().unwrap()).with_extension("o");
        let ok = Command::new(&cc)
            .args(["-std=c99", "-O2", "-fPIC", "-c", "-Ivendor/erfa"])
            .arg(src).arg("-o").arg(&obj)
            .status().expect("C compiler not found (set CC)").success();
        assert!(ok, "compiling {:?} failed", src);
        objects.push(obj);
        println!("cargo:rerun-if-changed={}", src.display());
    }
    let lib = out.join("liberfa_subset.a");
    let _ = fs::remove_file(&lib);
    let ok = Command::new(&ar).arg("crs").arg(&lib).args(&objects)
        .status().expect("ar not found (set AR)").success();
    assert!(ok, "ar failed");
    println!("cargo:rustc-link-search=native={}", out.display());
    println!("cargo:rustc-link-lib=static=erfa_subset");
    println!("cargo:rerun-if-changed=vendor/erfa/erfa.h");
    println!("cargo:rerun-if-changed=build.rs");
}
