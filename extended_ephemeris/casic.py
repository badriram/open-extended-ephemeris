# SPDX-License-Identifier: MIT
"""CASIC MSG-GPSEPH frames, for AT6558-based receivers (Quectel L76K, ATGM336H and others).

The 72-byte payload carries the LNAV ephemeris at its broadcast resolution; its first word,
"reserved" in CASIC's specification, is the sum of the payload's other 32-bit words. Verified
byte for byte against CASIC's own AGNSS file, and on an L76K, by
https://github.com/badriram/l76k-at6558-agnss, which also has a C encoder and casic_eph_select()
to pick each satellite's frame for the current time from a file of many."""
import struct

PI = 3.1415926535898                          # IS-GPS-200's pi, which semicircles are defined with


def checksum(cls, mid, payload):
    ck = (mid << 24) + (cls << 16) + len(payload)
    for i in range(0, len(payload) - len(payload) % 4, 4):
        ck += struct.unpack_from("<I", payload, i)[0]
    return ck & 0xFFFFFFFF


def frame(cls, mid, payload):
    return struct.pack("<BBHBB", 0xBA, 0xCE, len(payload), cls, mid) + payload + \
        struct.pack("<I", checksum(cls, mid, payload))


def split_frames(buf):
    """Check every frame's checksum; return {(class, id): count}."""
    j, kinds = 0, {}
    while j + 6 <= len(buf):
        ln = buf[j + 2] | buf[j + 3] << 8
        f = buf[j:j + 10 + ln]
        if f[:2] != b"\xba\xce" or struct.unpack("<I", f[-4:])[0] != checksum(f[4], f[5], f[6:-4]):
            raise ValueError(f"bad frame at offset {j}")
        kinds[(f[4], f[5])] = kinds.get((f[4], f[5]), 0) + 1
        j += len(f)
    return kinds


def ura_index(metres):
    """IS-GPS-200 URA index for an accuracy in metres."""
    for n, top in enumerate((2.4, 3.4, 4.85, 6.85, 9.65, 13.65, 24, 48, 96, 192, 384, 768, 1536, 3072, 6144)):
        if metres <= top:
            return n
    return 15


def _q(v, lsb, fmt):
    r = round(v / lsb)
    struct.pack("<" + fmt, r)                 # raises struct.error if out of range
    return r


def gpseph_payload(r, sv):
    """72-byte MSG-GPSEPH payload from a record with RINEX-style keys (angles in radians)."""
    s = 1 / PI
    body = struct.pack(
        "<IIIiiiiihhhhhhhhHHIihbbHBBBBH",
        0,
        _q(r["sqrtA"], 2**-19, "I"), _q(r["e"], 2**-33, "I"),
        _q(r["omega"] * s, 2**-31, "i"), _q(r["M0"] * s, 2**-31, "i"), _q(r["i0"] * s, 2**-31, "i"),
        _q(r["Omega0"] * s, 2**-31, "i"), _q(r["OmegaDot"] * s, 2**-43, "i"),
        _q(r["deltaN"] * s, 2**-43, "h"), _q(r["IDOT"] * s, 2**-43, "h"),
        _q(r["Cuc"], 2**-29, "h"), _q(r["Cus"], 2**-29, "h"), _q(r["Crc"], 2**-5, "h"),
        _q(r["Crs"], 2**-5, "h"), _q(r["Cic"], 2**-29, "h"), _q(r["Cis"], 2**-29, "h"),
        _q(r["toe"], 16, "H"), int(r["week"]) % 1024, _q(r["toc"], 16, "I"),
        _q(r["af0"], 2**-31, "i"), _q(r["af1"], 2**-43, "h"), _q(r["af2"], 2**-55, "b"),
        _q(r["TGD"], 2**-31, "b"), int(r["IODC"]), ura_index(r["svacc"]), int(r["health"]),
        sv, 3 if int(r["health"]) == 0 else 1, 0)
    word_sum = sum(struct.unpack_from("<17I", body, 4)) & 0xFFFFFFFF
    return struct.pack("<I", word_sum) + body[4:]


def gpseph_frame(r, sv):
    return frame(0x08, 0x07, gpseph_payload(r, sv))
