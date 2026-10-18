# SPDX-License-Identifier: MIT
"""A GPS orbit model for predicting days ahead.

Forces: Earth's gravity field to degree and order 4 (EGM2008), the Sun and the Moon (ERFA), and
solar radiation pressure as five ECOM-like terms (constant along the Sun direction, the solar
panel axis and the third axis, plus a once-per-revolution term on the third), switched off in
Earth's shadow. The equations of motion are integrated in GCRS; Earth orientation comes from ERFA
(IAU 2006/2000A, CIO-based) with polar motion and UT1 from the IGS ERP file. Treating the
Earth's axis as fixed in space instead left ~20 m of error over two days.
"""
import math

import erfa
import numpy as np
from scipy.integrate import solve_ivp

GM_E, R_E = 3.986004418e14, 6378136.3
GM_SUN, GM_MOON = 1.32712440018e20, 4.9028e12
AU = 1.495978707e11
JD_GPS0 = 2444244.5               # Julian date of the GPS epoch, 1980-01-06
TT_GPS = 51.184                   # TT - GPS, s
UTC_GPS = -18.0                   # UTC - GPS since 2017, s

# EGM2008, fully normalised. Degree 6 was tried and changed the fit by under 0.1 m.
_CBAR = {(2, 0): -0.484165143790815e-3, (2, 1): -0.206615509074176e-9, (2, 2): 0.243938357328313e-5,
         (3, 0): 0.957161207093473e-6, (3, 1): 0.203046201047864e-5, (3, 2): 0.904787894809528e-6,
         (3, 3): 0.721321757121568e-6, (4, 0): 0.539965866638991e-6, (4, 1): -0.536157389388867e-6,
         (4, 2): 0.350501623962649e-6, (4, 3): 0.990856766672321e-6, (4, 4): -0.188519633023033e-6}
_SBAR = {(2, 1): 0.138441389137979e-8, (2, 2): -0.140027370385934e-5, (3, 1): 0.248200415856872e-6,
         (3, 2): -0.619005475177618e-6, (3, 3): 0.141434926192941e-5, (4, 1): -0.473567346518086e-6,
         (4, 2): 0.662480026275829e-6, (4, 3): -0.200956723567452e-6, (4, 4): 0.308803882149194e-6}
NMAX = 4

# A priori radiation-pressure terms, nm/s^2: D0 (along the Sun direction) is left free; the others
# are usually a few nm/s^2 for GPS and are held near zero unless the data insists.
SRP_PRIOR_MEAN = np.array([-100.0, 0.0, 0.0, 0.0, 0.0])
SRP_PRIOR_SIGMA = np.array([1e6, 2.0, 5.0, 5.0, 5.0])
_C = np.zeros((NMAX + 1, NMAX + 1))
_S = np.zeros((NMAX + 1, NMAX + 1))
_C[0, 0] = 1.0
for (_n, _m), _v in _CBAR.items():
    _f = math.sqrt((2 - (_m == 0)) * (2 * _n + 1) * math.factorial(_n - _m) / math.factorial(_n + _m))
    _C[_n, _m] = _v * _f
    _S[_n, _m] = _SBAR.get((_n, _m), 0.0) * _f


def gravity(r):
    """Geopotential acceleration in the Earth-fixed frame (Montenbruck & Gill, eq. 3.33)."""
    x, y, z = r
    r2 = x * x + y * y + z * z
    rho = R_E * R_E / r2
    x0, y0, z0 = R_E * x / r2, R_E * y / r2, R_E * z / r2
    N = NMAX + 2
    V = np.zeros((N, N))
    W = np.zeros((N, N))
    V[0, 0] = R_E / math.sqrt(r2)
    V[1, 0] = z0 * V[0, 0]
    for n in range(2, N):
        V[n, 0] = ((2 * n - 1) * z0 * V[n - 1, 0] - (n - 1) * rho * V[n - 2, 0]) / n
    for m in range(1, N):
        V[m, m] = (2 * m - 1) * (x0 * V[m - 1, m - 1] - y0 * W[m - 1, m - 1])
        W[m, m] = (2 * m - 1) * (x0 * W[m - 1, m - 1] + y0 * V[m - 1, m - 1])
        if m + 1 < N:
            V[m + 1, m] = (2 * m + 1) * z0 * V[m, m]
            W[m + 1, m] = (2 * m + 1) * z0 * W[m, m]
        for n in range(m + 2, N):
            V[n, m] = ((2 * n - 1) * z0 * V[n - 1, m] - (n + m - 1) * rho * V[n - 2, m]) / (n - m)
            W[n, m] = ((2 * n - 1) * z0 * W[n - 1, m] - (n + m - 1) * rho * W[n - 2, m]) / (n - m)
    ax = ay = az = 0.0
    for n in range(NMAX + 1):
        for m in range(n + 1):
            C, S = _C[n, m], _S[n, m]
            if m == 0:
                ax -= C * V[n + 1, 1]
                ay -= C * W[n + 1, 1]
                az -= (n + 1) * C * V[n + 1, 0]
            else:
                f = 0.5 * (n - m + 1) * (n - m + 2)
                ax += 0.5 * (-C * V[n + 1, m + 1] - S * W[n + 1, m + 1]) + f * (C * V[n + 1, m - 1] + S * W[n + 1, m - 1])
                ay += 0.5 * (-C * W[n + 1, m + 1] + S * V[n + 1, m + 1]) + f * (-C * W[n + 1, m - 1] + S * V[n + 1, m - 1])
                az += (n - m + 1) * (-C * V[n + 1, m] - S * W[n + 1, m])
    return GM_E / (R_E * R_E) * np.array([ax, ay, az])


def gravity_many(R):
    """gravity() for k positions at once: R is (k, 3), the result (k, 3)."""
    R = np.atleast_2d(R)
    x, y, z = R[:, 0], R[:, 1], R[:, 2]
    r2 = x * x + y * y + z * z
    rho = R_E * R_E / r2
    x0, y0, z0 = R_E * x / r2, R_E * y / r2, R_E * z / r2
    N, k = NMAX + 2, len(x)
    V = np.zeros((N, N, k))
    W = np.zeros((N, N, k))
    V[0, 0] = R_E / np.sqrt(r2)
    V[1, 0] = z0 * V[0, 0]
    for n in range(2, N):
        V[n, 0] = ((2 * n - 1) * z0 * V[n - 1, 0] - (n - 1) * rho * V[n - 2, 0]) / n
    for m in range(1, N):
        V[m, m] = (2 * m - 1) * (x0 * V[m - 1, m - 1] - y0 * W[m - 1, m - 1])
        W[m, m] = (2 * m - 1) * (x0 * W[m - 1, m - 1] + y0 * V[m - 1, m - 1])
        if m + 1 < N:
            V[m + 1, m] = (2 * m + 1) * z0 * V[m, m]
            W[m + 1, m] = (2 * m + 1) * z0 * W[m, m]
        for n in range(m + 2, N):
            V[n, m] = ((2 * n - 1) * z0 * V[n - 1, m] - (n + m - 1) * rho * V[n - 2, m]) / (n - m)
            W[n, m] = ((2 * n - 1) * z0 * W[n - 1, m] - (n + m - 1) * rho * W[n - 2, m]) / (n - m)
    ax, ay, az = np.zeros(k), np.zeros(k), np.zeros(k)
    for n in range(NMAX + 1):
        for m in range(n + 1):
            C, S = _C[n, m], _S[n, m]
            if m == 0:
                ax -= C * V[n + 1, 1]
                ay -= C * W[n + 1, 1]
                az -= (n + 1) * C * V[n + 1, 0]
            else:
                f = 0.5 * (n - m + 1) * (n - m + 2)
                ax += 0.5 * (-C * V[n + 1, m + 1] - S * W[n + 1, m + 1]) + f * (C * V[n + 1, m - 1] + S * W[n + 1, m - 1])
                ay += 0.5 * (-C * W[n + 1, m + 1] + S * V[n + 1, m + 1]) + f * (-C * W[n + 1, m - 1] + S * V[n + 1, m - 1])
                az += (n - m + 1) * (-C * V[n + 1, m] - S * W[n + 1, m])
    return GM_E / (R_E * R_E) * np.stack([ax, ay, az], axis=1)


class EarthOrientation:
    """GCRS <-> ITRS, and the Sun and Moon, from ERFA and IGS ERP rows (gps seconds, xp rad, yp rad,
    UT1-UTC s). Everything that changes slowly is tabulated every 30 minutes and interpolated:
    calling ERFA at every integration step was most of the run time."""

    STEP = 1800.0

    def __init__(self, erp, t_start, t_end):
        from scipy.interpolate import CubicSpline
        ts = [row[0] for row in erp]
        self.grid = np.arange(t_start - 3 * 3600, t_end + 86400, self.STEP)
        self._dut1 = np.interp(self.grid, ts, [row[3] for row in erp])
        mats, sun, moon = [], [], []
        for t, xp, yp in zip(self.grid, np.interp(self.grid, ts, [row[1] for row in erp]),
                             np.interp(self.grid, ts, [row[2] for row in erp])):
            d1, d2 = JD_GPS0, (t + TT_GPS) / 86400.0
            mats.append(erfa.pom00(xp, yp, erfa.sp00(d1, d2)).ravel())    # polar motion ...
            # CIO-based celestial-to-intermediate matrix, which pairs with the Earth rotation angle.
            # (The equinox-based pnm06a pairs with sidereal time instead: mixing it with the ERA
            # turned the Earth by ~0.5 degrees.)
            mats.append(erfa.c2i06a(d1, d2).ravel())
            pvh, _ = erfa.epv00(d1, d2)
            sun.append(-pvh[0] * AU)
            moon.append(erfa.moon98(d1, d2)[0] * AU)
        mats = np.array(mats).reshape(len(self.grid), 2, 9)
        self._pom, self._c2i = mats[:, 0], mats[:, 1]
        self._sun, self._moon = CubicSpline(self.grid, np.array(sun)), CubicSpline(self.grid, np.array(moon))

    def _lin(self, table, t):
        k = min(max(int((t - self.grid[0]) // self.STEP), 0), len(self.grid) - 2)
        f = (t - self.grid[k]) / self.STEP
        return table[k] * (1 - f) + table[k + 1] * f

    def c2t(self, t):
        dut1 = self._lin(self._dut1, t)
        du = JD_GPS0 - 2451545.0 + (t + UTC_GPS + dut1) / 86400.0          # UT1 days from J2000
        era = 2 * math.pi * ((0.7790572732640 + 0.00273781191135448 * du + du) % 1.0)
        c, s = math.cos(era), math.sin(era)
        rz = np.array([[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]])
        return self._lin(self._pom, t).reshape(3, 3) @ rz @ self._lin(self._c2i, t).reshape(3, 3)

    def bodies(self, t):
        return self._sun(t), self._moon(t)


def read_erp(text):
    """IGS ERP file -> rows (gps seconds, xp rad, yp rad, UT1-UTC s), extended ten days ahead with
    UT1 drifting by the last length of day and the pole held."""
    arc = math.pi / 180 / 3600 * 1e-6
    rows = []
    for line in text.splitlines():
        f = line.split()
        if f and f[0][0].isdigit() and len(f) >= 5:
            t = (float(f[0]) - 44244.0) * 86400 - UTC_GPS
            rows.append((t, float(f[1]) * arc, float(f[2]) * arc, float(f[3]) * 1e-7, float(f[4]) * 1e-7))
    if not rows:
        raise ValueError("no ERP rows")
    t, xp, yp, ut, lod = rows[-1]
    rows += [(t + d * 86400, xp, yp, ut - lod * d, lod) for d in range(1, 11)]
    return [row[:4] for row in rows]


def _lagrange(ts, ys, t, order=10):
    i = int(np.searchsorted(ts, t))
    lo = max(0, min(len(ts) - order, i - order // 2))
    xs = ts[lo:lo + order]
    out = 0
    for j in range(order):
        w = 1.0
        for m in range(order):
            if m != j:
                w *= (t - xs[m]) / (xs[j] - xs[m])
        out = out + w * ys[lo + j]
    return out


R_SUN = 6.957e8
# Integration: 8th-order Dormand-Prince with a fixed 120 s step, which agrees with an adaptive
# 1e-14 run to 12 mm over 48 h. Adaptive steps were noisy by centimetres to metres (at 1e-11 the
# orbit was ~2.7 m off after 48 h, by amounts that jumped with the steps chosen).
STEP_S = 120.0


def sunlit(r, sun):
    """Fraction of sunlight reaching a satellite at r: 1 in sunlight, 0 in Earth's umbra, and a smooth
    ramp across the penumbra (~190 km wide at GPS distance). Switching radiation pressure off
    abruptly made the integration depend on where its steps fell by metres, and the fit with it."""
    sn = np.linalg.norm(sun)
    ss = sun / sn
    proj = r @ ss
    if proj >= 0:
        return 1.0
    rho = np.linalg.norm(r - proj * ss)                # distance from the shadow's axis
    w = 2 * (-proj) * R_SUN / sn                       # penumbra width at this distance behind Earth
    x = (rho - (R_E - w / 2)) / w
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0
    return x * x * (3 - 2 * x)                         # smoothstep: continuous, with a continuous slope


class Orbit:
    def __init__(self, eo, t0):
        self.eo, self.t0 = eo, t0

    def _bodies(self, t):
        return self.eo.bodies(t)

    def _accel_many(self, t, R, V, P):
        """Accelerations for k trajectories at time t: R, V (k, 3) positions and velocities in GCRS,
        P (k, 5) radiation-pressure terms in nm/s^2."""
        m = self.eo.c2t(t)
        a = gravity_many(R @ m.T) @ m
        sun, moon = self._bodies(t)
        for s, gm in ((sun, GM_SUN), (moon, GM_MOON)):
            d = s - R
            a += gm * (d / np.linalg.norm(d, axis=1)[:, None] ** 3 - s / np.linalg.norm(s) ** 3)
        sn = np.linalg.norm(sun)
        ss = sun / sn
        lit = np.array([sunlit(r, sun) for r in R])
        if lit.any():
            eD = sun - R
            eD /= np.linalg.norm(eD, axis=1)[:, None]
            er = R / np.linalg.norm(R, axis=1)[:, None]
            eY = np.cross(eD, er)
            eY /= np.linalg.norm(eY, axis=1)[:, None]
            eB = np.cross(eD, eY)
            h = np.cross(R, V)
            h /= np.linalg.norm(h, axis=1)[:, None]
            sp = ss - (h @ ss)[:, None] * h
            sp /= np.linalg.norm(sp, axis=1)[:, None]
            u = np.arctan2(np.einsum("ij,ij->i", np.cross(sp, er), h), np.einsum("ij,ij->i", sp, er))
            Bt = P[:, 2] + P[:, 3] * np.cos(u) + P[:, 4] * np.sin(u)
            a += 1e-9 * lit[:, None] * (P[:, 0, None] * eD + P[:, 1, None] * eY + Bt[:, None] * eB)
        return a

    def _rhs_many(self, t, y, P):
        Y = y.reshape(-1, 6)
        A = self._accel_many(t, Y[:, :3], Y[:, 3:], P)
        return np.concatenate([Y[:, 3:], A], axis=1).ravel()

    def propagate_many(self, states, srps, times):
        """Positions (k, len(times), 3) for k trajectories, integrated together with fixed steps.
        Fixed steps make the result a smooth function of the inputs: with adaptive steps a
        0.01 nm/s^2 change in radiation pressure "moved" the orbit by half a metre of noise, and
        the fit could not converge on that."""
        states, srps = np.atleast_2d(states), np.atleast_2d(srps)
        sol = solve_ivp(self._rhs_many, (self.t0, float(max(times))), states.ravel(), method="DOP853",
                        t_eval=times, first_step=STEP_S, max_step=STEP_S, rtol=1e3, atol=1e12,
                        args=(srps,))
        return sol.y.reshape(len(states), 6, -1)[:, :3, :].transpose(0, 2, 1)

    def propagate(self, state, srp, times):
        return self.propagate_many(state, srp, times)[0]

    def state_at(self, state, srp, t):
        """The full state (position and velocity) at time t, from state at t0."""
        sol = solve_ivp(self._rhs_many, (self.t0, float(t)), np.asarray(state, float), method="DOP853",
                        first_step=STEP_S, max_step=STEP_S, rtol=1e3, atol=1e12,
                        args=(np.atleast_2d(srp),))
        return sol.y[:, -1]

    # finite-difference steps for the Jacobian: metres, metres per second, nm/s^2
    _STEPS = np.array([10.0, 10.0, 10.0, 1e-2, 1e-2, 1e-2, 10.0, 10.0, 10.0, 10.0, 10.0])

    def fit(self, times, positions_gcrs, start=None, sigma_obs=1.0, max_iter=25):
        """Initial state and radiation-pressure terms that best reproduce the positions: Levenberg-
        Marquardt on scaled unknowns, run until a step moves the orbit by under a millimetre, with a Jacobian by finite
        differences (exact enough now that the integration is clean, where cheaper analytic
        partials stopped metres short). The secondary radiation-pressure terms carry a priori
        constraints (SRP_PRIOR_*). Returns (state, srp, rms m, iterations, converged). start: an
        optional (state, srp) to begin from; it changes how quickly the fit gets there, not where."""
        ts, xs = np.asarray(times, float), np.asarray(positions_gcrs, float)
        if start is None:
            v0 = (_lagrange(ts, xs, ts[0] + 1.0) - _lagrange(ts, xs, ts[0] - 1.0)) / 2.0
            p = np.concatenate([xs[0], v0, SRP_PRIOR_MEAN])
        else:
            p = np.concatenate([start[0], start[1]])
        lam = np.concatenate([np.zeros(6), 1.0 / np.asarray(SRP_PRIOR_SIGMA) ** 2])
        prior = np.concatenate([np.zeros(6), SRP_PRIOR_MEAN])

        def cost(q, pos):
            return float(np.sum((pos - xs) ** 2)) / sigma_obs ** 2 + float(np.sum(lam * (q - prior) ** 2))

        base = self.propagate(p[:6], p[6:], ts)
        c = cost(p, base)
        mu, converged = 1e-3, False
        for it in range(1, max_iter + 1):
            Q = np.tile(p, (11, 1)) + np.diag(self._STEPS)
            many = self.propagate_many(Q[:, :6], Q[:, 6:], ts)
            J = np.array([((many[j] - base) / self._STEPS[j]).ravel() for j in range(11)]).T
            res = (base - xs).ravel()
            # Scale every unknown to unit effect: velocity partials are ~1e5 times position ones, and
            # unscaled the normal matrix (condition ~1e14) turned rounding into steps of decimetres.
            d = np.linalg.norm(J, axis=0)
            Js = J / d
            with np.errstate(all="ignore"):     # numpy 2.0 + Accelerate: spurious matmul warnings
                N = Js.T @ Js / sigma_obs ** 2 + np.diag(lam / d ** 2)
                g = Js.T @ res / sigma_obs ** 2 + lam * (p - prior) / d
            while True:                          # Levenberg-Marquardt: damp until the cost falls
                with np.errstate(all="ignore"):
                    dp = -np.linalg.solve(N + mu * np.diag(np.diag(N)), g) / d
                    moved = float(np.max(np.abs(J @ dp)))
                q = p + dp
                pos = self.propagate(q[:6], q[6:], ts)
                cq = cost(q, pos)
                if cq <= c or moved < 1e-5:
                    break
                mu *= 10
                if mu > 1e6:
                    break
            if cq <= c:
                p, base, c, mu = q, pos, cq, max(mu / 10, 1e-6)
            if moved < 1e-4:                     # the step moves the orbit by under 0.1 mm
                converged = True
                break
        return p[:6], p[6:], math.sqrt(float(np.mean((base - xs) ** 2))), it, converged
