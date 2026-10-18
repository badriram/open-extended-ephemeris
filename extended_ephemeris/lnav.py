# SPDX-License-Identifier: MIT
"""GPS LNAV broadcast ephemeris (IS-GPS-200): evaluate one, and fit one to a trajectory.

Positions are ECEF metres; times are continuous GPS seconds (since 1980-01-06 00:00 GPS)."""
import math

import numpy as np
from scipy.optimize import least_squares

GM = 3.986005e14                 # IS-GPS-200's values, which receivers use
OMEGA_E = 7.2921151467e-5
WEEK = 604800.0
KEYS = ["sqrtA", "e", "M0", "omega", "i0", "Omega0", "OmegaDot", "deltaN", "IDOT",
        "Cuc", "Cus", "Crc", "Crs", "Cic", "Cis"]


def sat_pos(p, toe, t):
    """ECEF position at t from ephemeris parameters p (in KEYS order) with reference time toe."""
    sqrtA, e, M0, w, i0, W0, Wd, dn, idot, cuc, cus, crc, crs, cic, cis = p
    A = sqrtA * sqrtA
    tk = t - toe
    M = M0 + (math.sqrt(GM / A ** 3) + dn) * tk
    E = M
    for _ in range(12):
        E -= (E - e * math.sin(E) - M) / (1 - e * math.cos(E))
    phi = math.atan2(math.sqrt(1 - e * e) * math.sin(E), math.cos(E) - e) + w
    s2, c2 = math.sin(2 * phi), math.cos(2 * phi)
    u = phi + cus * s2 + cuc * c2
    r = A * (1 - e * math.cos(E)) + crs * s2 + crc * c2
    i = i0 + idot * tk + cis * s2 + cic * c2
    x, y = r * math.cos(u), r * math.sin(u)
    W = W0 + (Wd - OMEGA_E) * tk - OMEGA_E * (toe % WEEK)
    cW, sW, ci = math.cos(W), math.sin(W), math.cos(i)
    return np.array([x * cW - y * ci * sW, x * sW + y * ci * cW, y * math.sin(i)])


def _elements(r_ecef, v_ecef, toe):
    """Osculating elements, in broadcast form, from an ECEF state at toe: the fit's starting point."""
    r = np.asarray(r_ecef, float)
    v = np.asarray(v_ecef, float) + np.cross([0, 0, OMEGA_E], r)
    h = np.cross(r, v)
    hn = np.linalg.norm(h)
    i = math.acos(h[2] / hn)
    node = np.cross([0, 0, 1], h)
    nhat = node / np.linalg.norm(node)
    ev = np.cross(v, h) / GM - r / np.linalg.norm(r)
    e = float(np.linalg.norm(ev))
    a = 1 / (2 / np.linalg.norm(r) - v @ v / GM)
    w = math.atan2(np.cross(nhat, ev) @ h / hn, nhat @ ev)
    u = math.atan2(np.cross(nhat, r) @ h / hn, nhat @ r)
    nu = u - w
    E = math.atan2(math.sqrt(1 - e * e) * math.sin(nu), e + math.cos(nu))
    W0 = (math.atan2(node[1], node[0]) + OMEGA_E * (toe % WEEK) + math.pi) % (2 * math.pi) - math.pi
    return np.array([math.sqrt(a), e, E - e * math.sin(E), w, i, W0, -8.0e-9] + [0.0] * 8)


# The fit runs in non-singular variables (e cos w, e sin w, M0 + w), which stay well defined for
# GPS's near-circular orbits, scaled so that a unit step moves the satellite by roughly a metre.
_SCALE = np.array([1e-4, 1e-7, 1e-8, 1e-7, 1e-8, 1e-8, 1e-11, 1e-11, 1e-11,
                   1e-7, 1e-7, 1, 1, 1e-7, 1e-7])


def _to_ns(p):
    q = p.copy()
    q[1], q[3], q[2] = p[1] * math.cos(p[3]), p[1] * math.sin(p[3]), p[2] + p[3]
    return q


def _from_ns(q):
    p = q.copy()
    p[1], p[3] = math.hypot(q[1], q[3]), math.atan2(q[3], q[1])
    p[2] = (q[2] - p[3] + math.pi) % (2 * math.pi) - math.pi
    p[5] = (p[5] + math.pi) % (2 * math.pi) - math.pi
    return p


def fit(times, positions, toe):
    """Fit ephemeris parameters to positions (ECEF, m) at times spanning toe +/- 2 h.
    Returns (parameters in KEYS order, rms residual in metres)."""
    ts, xs = np.asarray(times, float), np.asarray(positions, float)
    k = min(max(int(np.argmin(abs(ts - toe))), 1), len(ts) - 2)
    p0 = _elements(xs[k], (xs[k + 1] - xs[k - 1]) / (ts[k + 1] - ts[k - 1]), ts[k])
    p0[2] += math.sqrt(GM / p0[0] ** 6) * (toe - ts[k])
    q0 = _to_ns(p0)

    def res(d):
        p = _from_ns(q0 + d * _SCALE)
        return np.concatenate([sat_pos(p, toe, t) - x for t, x in zip(ts, xs)])

    def jac(d, h=1e-2):
        # Central differences with ~1 cm steps. The solver's default steps are so small that the
        # change is lost in rounding at GPS distances and the fit wanders off.
        cols = []
        for j in range(15):
            step = np.zeros(15)
            step[j] = h
            cols.append((res(d + step) - res(d - step)) / (2 * h))
        return np.array(cols).T

    sol = least_squares(res, np.zeros(15), jac=jac, method="lm", xtol=1e-12, ftol=1e-12, max_nfev=200)
    return _from_ns(q0 + sol.x * _SCALE), math.sqrt(float(np.mean(sol.fun ** 2)))
