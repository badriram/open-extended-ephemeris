# Predicted GPS and Galileo orbit data (extended ephemeris)

A week of GPS and Galileo broadcast-style ephemerides predicted from the newest CODE 5-day orbit by
https://github.com/badriram/open-extended-ephemeris. These are predictions, not the
satellites' own broadcasts. gps-7d.rnx is RINEX 3.04; gps-7d-ubx.bin is u-blox UBX-MGA-GPS-EPH
and UBX-MGA-GAL-EPH frames (M8/M9/M10); gps-7d-casic.bin is GPS-only CASIC MSG-GPSEPH frames for AT6558 receivers (Quectel
L76K). "extended-ephemeris select" picks the frames for now. gps-7d.json says when and from what it was made,
with a self-check against the satellites' broadcasts. This branch holds only the latest.

Orbits: Center for Orbit Determination in Europe (CODE), University of Bern, DOI 10.48620/97774;
International GNSS Service (IGS), via BKG. Broadcast clocks: IGS, via BKG. MIT licence (see LICENSE).
