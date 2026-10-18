# ERFA subset

37 source files and two headers from ERFA 2.0.1 (https://github.com/liberfa/erfa, release
tarball sha256 45111f0addf6ef69345aaf093bd054705b41ac2db7541dd3f0e105362c8e20e4), unmodified:
the routines behind eraC2i06a, eraEra00, eraPom00, eraSp00, eraEpv00 and eraMoon98 and everything
they call. `build.rs` compiles them with the system C compiler. ERFA is derived from the IAU's SOFA
library; its BSD-3-Clause licence is in LICENSE here and applies to these files.
