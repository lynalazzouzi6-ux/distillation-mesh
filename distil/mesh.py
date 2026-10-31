"""
Simulateur MESH rigoureux d'une colonne de distillation binaire (étages théoriques, P constante,
condenseur total, rebouilleur partiel = dernier étage, une alimentation).

Étages numérotés 1..N de la tête vers le pied (étage N = rebouilleur). Dans le code (indices 0..N-1).
Spécifications : distillat D (kmol/h) et taux de reflux R = L0/D.
Inconnues (4N+1) : x_j, T_j, L_j (j=1..N), V_j (j=2..N), Q_R (chaleur du rebouilleur), T_0 (T du reflux).
Équations (4N+1), pour chaque étage j :
  M (bilan matière total)         L_{j-1} + V_{j+1} + F_j - L_j - V_j = 0
  M (bilan composant léger)       L_{j-1} x_{j-1} + V_{j+1} y_{j+1} + F_j zF - L_j x_j - V_j y_j = 0
  E+S (équilibre + sommation)     K_A x_j + K_B (1 - x_j) = 1        avec y_j = K_A x_j
  H (bilan enthalpique)           L_{j-1} h_{j-1} + V_{j+1} H_{j+1} + F_j hF - L_j h_j - V_j H_j + Q_j = 0
+ équation de bulle du reflux (condenseur total : x_0 = y_1, liquide saturé à T_0).
Résolution par Newton amorti ; jacobien par différences finies "colorées" (structure tridiagonale
par blocs : 14 évaluations du résidu par itération, quel que soit N).
Efficacité : le MESH utilise des étages théoriques ; on passe aux plateaux réels par N_réel = N / E_o.
"""
from dataclasses import dataclass, field
import numpy as np


@dataclass
class SolutionMESH:
    N: int; NF: int; F: float; zF: float; q: float; D: float; W: float; R: float
    x: np.ndarray; y: np.ndarray; T: np.ndarray; L: np.ndarray; V: np.ndarray
    T0: float; Qc: float; Qr: float          # MJ/h
    xD: float; xW: float
    convergee: bool; iterations: int; residu: float
    u: np.ndarray = field(repr=False, default=None)
    bilans: dict = field(default_factory=dict)


def _hF(sys, zF, q):
    return sys.enthalpie_alimentation(zF, q)


class _Probleme:
    def __init__(self, sys, N, NF, F, zF, q, D, R):
        self.s, self.N, self.NF = sys, N, NF
        self.F, self.zF, self.q, self.D, self.R = F, zF, q, D, R
        self.W = F - D
        self.hF = _hF(sys, zF, q)
        self.V1 = (R + 1.0) * D
        self.sc_en = F * 40.0                      # échelle des équations d'énergie

    # --- dépaquetage
    def unpack(self, u):
        N = self.N
        x, T, L = u[0:N], u[N:2 * N], u[2 * N:3 * N]
        V = np.concatenate(([self.V1], u[3 * N:4 * N - 1]))
        return x, T, L, V, u[4 * N - 1], u[4 * N]

    def residu(self, u):
        s, N, NF = self.s, self.N, self.NF
        x, T, L, V, QR, T0 = self.unpack(u)
        K1, K2 = s.K(x, T)
        y = K1 * x
        h = s.h_liq(x, T)
        H = s.H_vap(y, T)
        # reflux : liquide saturé de composition y[0] à T0
        k1_0, k2_0 = s.K(y[0], T0)
        h0 = s.h_liq(y[0], T0)
        Lin = np.concatenate(([self.R * self.D], L[:-1]))
        xin = np.concatenate(([y[0]], x[:-1]))
        hin = np.concatenate(([h0], h[:-1]))
        Vin = np.concatenate((V[1:], [0.0]))
        yin = np.concatenate((y[1:], [0.0]))
        Hin = np.concatenate((H[1:], [0.0]))
        Ff = np.zeros(N); Ff[NF - 1] = self.F
        Q = np.zeros(N); Q[N - 1] = QR
        f_tot = (Lin + Vin + Ff - L - V) / self.F
        f_cmp = (Lin * xin + Vin * yin + Ff * self.zF - L * x - V * y) / self.F
        f_bub = K1 * x + K2 * (1.0 - x) - 1.0
        f_en = (Lin * hin + Vin * Hin + Ff * self.hF - L * h - V * H + Q) / self.sc_en
        f_T0 = k1_0 * y[0] + k2_0 * (1.0 - y[0]) - 1.0
        return np.concatenate((f_tot, f_cmp, f_bub, f_en, [f_T0]))

    # --- jacobien par différences finies colorées
    def jacobien(self, u, f0):
        N = self.N
        n = 4 * N + 1
        J = np.zeros((n, n))
        rows_stage = lambda k: [k, N + k, 2 * N + k, 3 * N + k]     # lignes des 4 équations de l'étage k

        def col(typ, j):
            if typ == 0: return j
            if typ == 1: return N + j
            if typ == 2: return 2 * N + j
            return 3 * N + j - 1                                    # V_j, j >= 1

        for typ in range(4):
            for c in range(3):
                idx = [j for j in range(N) if j % 3 == c and not (typ == 3 and j == 0)]
                if not idx: continue
                up = u.copy(); hs = {}
                for j in idx:
                    k = col(typ, j)
                    hs[j] = 1e-6 * max(1.0, abs(u[k])); up[k] += hs[j]
                df = self.residu(up) - f0
                for j in idx:
                    k = col(typ, j)
                    for st in (j - 1, j, j + 1):
                        if 0 <= st < N:
                            for r in rows_stage(st):
                                J[r, k] = df[r] / hs[j]
                    if j == 0 and typ in (0, 1):                    # ligne de l'équation de T0
                        J[4 * N, k] = df[4 * N] / hs[j]
        for k in (4 * N - 1, 4 * N):                                # Q_R et T0 : colonnes pleines
            up = u.copy(); hh = 1e-6 * max(1.0, abs(u[k])); up[k] += hh
            J[:, k] = (self.residu(up) - f0) / hh
        return J


def estimation_initiale(sys, N, NF, F, zF, q, D, R, xD_guess, xW_guess):
    """Profil initial : x linéaire par morceaux, débits à la McCabe-Thiele (CMO)."""
    W = F - D
    j = np.arange(N)
    x = np.where(j <= NF - 1,
                 xD_guess + (zF - xD_guess) * j / max(NF - 1, 1),
                 zF + (xW_guess - zF) * (j - (NF - 1)) / max(N - 1 - (NF - 1), 1))
    x = np.clip(x, 1e-4, min(sys.x_max, 0.9999) - 1e-4)
    T = np.array([sys.T_bulle(xx) for xx in x])
    Vr = (R + 1) * D
    Vs = Vr - (1 - q) * F
    L = np.where(j <= NF - 2, R * D, R * D + q * F).astype(float)
    L[N - 1] = W
    V = np.where(j <= NF - 1, Vr, Vs).astype(float)
    y0 = sys.y_eq(x[0])
    T0 = sys.T_bulle(y0)
    # chaleur du rebouilleur estimée par le bilan du dernier étage
    y = np.array([sys.y_eq(xx) for xx in x])
    h = sys.h_liq(x, T); H = sys.H_vap(y, T)
    QR = V[N - 1] * H[N - 1] + W * h[N - 1] - L[N - 2] * h[N - 2] if N > 1 else 0.0
    return np.concatenate((x, T, L, V[1:], [QR, T0]))


def resoudre(sys, N, NF, F, zF, q, D, R, u0=None, xD_guess=None, xW_guess=None,
             tol=1e-9, max_iter=80):
    """Résout le MESH. Retourne SolutionMESH (convergee=False si échec)."""
    if not (2 <= NF <= N - 1) and N > 2:
        NF = min(max(NF, 2), N - 1)
    P = _Probleme(sys, N, NF, F, zF, q, D, R)
    if u0 is None:
        xDg = xD_guess if xD_guess else min(0.9 * sys.x_max, 0.95)
        xWg = xW_guess if xW_guess else 0.02
        u = estimation_initiale(sys, N, NF, F, zF, q, D, R, xDg, xWg)
    else:
        u = u0.copy()
    f = P.residu(u)
    nrm = np.max(np.abs(f))
    it = 0
    for it in range(1, max_iter + 1):
        if not np.isfinite(nrm): break
        if nrm < tol: it -= 1; break
        J = P.jacobien(u, f)
        try:
            du = np.linalg.solve(J, -f)
        except np.linalg.LinAlgError:
            break
        # limitation de pas : x, T, flux
        lam = 1.0
        mx = np.max(np.abs(du[0:N]))
        if mx > 0.25: lam = 0.25 / mx
        mt = np.max(np.abs(du[N:2 * N]))
        if mt > 25.0: lam = min(lam, 25.0 / mt)
        accepte = False
        for _ in range(12):
            un = u + lam * du
            un[0:N] = np.clip(un[0:N], 1e-8, 1 - 1e-8)
            try:
                fn = P.residu(un)
                nn = np.max(np.abs(fn))
            except FloatingPointError:
                nn = np.inf
            if np.isfinite(nn) and nn < nrm * (1 - 1e-4 * lam) + 1e-14:
                accepte = True; break
            lam *= 0.5
        if not accepte: break
        u, f, nrm = un, fn, nn
    conv = bool(np.isfinite(nrm) and nrm < 1e-7)

    x, T, L, V, QR, T0 = P.unpack(u)
    K1, K2 = sys.K(x, T)
    y = K1 * x
    h = sys.h_liq(x, T); H = sys.H_vap(y, T)
    h0 = float(sys.h_liq(y[0], T0))
    Qc = float(V[0] * H[0] - (R * D + D) * h0)
    # bilans globaux (doivent se fermer)
    bilan = dict(
        matiere=abs(F - D - L[-1]),
        composant=abs(F * zF - D * y[0] - L[-1] * x[-1]),
        energie=abs(F * P.hF + QR - Qc - D * h0 - L[-1] * h[-1]) / (F * 40.0),
    )
    return SolutionMESH(N=N, NF=NF, F=F, zF=zF, q=q, D=D, W=F - D, R=R, x=x, y=y, T=T, L=L, V=V,
                        T0=float(T0), Qc=Qc, Qr=float(QR), xD=float(y[0]), xW=float(x[-1]),
                        convergee=conv, iterations=it, residu=float(nrm), u=u.copy(), bilans=bilan)
