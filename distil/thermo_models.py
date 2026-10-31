"""
Thermodynamique d'un binaire à P constante :
  - coefficients d'activité : NRTL ou UNIQUAC (vectorisés numpy)
  - équilibre liquide-vapeur : loi de Raoult modifiée  y_i P = x_i γ_i Psat_i
  - enthalpies : h(x,T) liquide (avec enthalpie d'excès) et H(y,T) vapeur
Convention : composant A = le plus volatil (léger), B = le plus lourd.
x, y = fractions molaires de A.  T en °C.  Enthalpies en kJ/mol.
"""
import numpy as np
from scipy.optimize import brentq
from scipy.interpolate import PchipInterpolator
from .components import COMPOSANTS

R_GAZ = 8.314462e-3   # kJ/mol/K
T_REF = 0.0           # °C : référence des enthalpies (liquides purs à 0 °C)


class Systeme:
    def __init__(self, cle, A, B, modele, params, P_mmHg=760.0, defauts=None,
                 donnees_exp=None, description=""):
        self.cle, self.modele, self.params = cle, modele, params
        self.cA, self.cB = COMPOSANTS[A], COMPOSANTS[B]
        self.nomA, self.nomB = self.cA["nom"], self.cB["nom"]
        self.P = P_mmHg
        self.defauts = defauts or {}
        self.donnees_exp = donnees_exp
        self.description = description
        self._construire_tables()

    # ------------------------------------------------------------------ pressions de vapeur
    @staticmethod
    def _psat(T, ant):
        A, B, C = ant
        return 10.0 ** (A - B / (T + C))

    def psat(self, T):
        return self._psat(T, self.cA["antoine"]), self._psat(T, self.cB["antoine"])

    # ------------------------------------------------------------------ coefficients d'activité
    def gammas(self, x, T):
        """Retourne (gamma_A, gamma_B). x, T : scalaires ou tableaux numpy. T en °C."""
        x = np.clip(np.asarray(x, dtype=float), 1e-12, 1 - 1e-12)
        TK = np.asarray(T, dtype=float) + 273.15
        x1, x2 = x, 1.0 - x
        if self.modele == "NRTL":
            b12, b21, al = self.params["b12"], self.params["b21"], self.params["alpha"]
            t12, t21 = b12 / TK, b21 / TK
            G12, G21 = np.exp(-al * t12), np.exp(-al * t21)
            ln1 = x2**2 * (t21 * (G21 / (x1 + x2 * G21))**2 + t12 * G12 / (x2 + x1 * G12)**2)
            ln2 = x1**2 * (t12 * (G12 / (x2 + x1 * G12))**2 + t21 * G21 / (x1 + x2 * G21)**2)
            return np.exp(ln1), np.exp(ln2)
        if self.modele == "UNIQUAC":
            r1, r2, q1, q2 = self.cA["r"], self.cB["r"], self.cA["q"], self.cB["q"]
            t12 = np.exp(self.params["b12"] / TK)
            t21 = np.exp(self.params["b21"] / TK)
            phi1 = r1 * x1 / (r1 * x1 + r2 * x2); phi2 = 1.0 - phi1
            th1 = q1 * x1 / (q1 * x1 + q2 * x2); th2 = 1.0 - th1
            l1 = 5.0 * (r1 - q1) - (r1 - 1.0); l2 = 5.0 * (r2 - q2) - (r2 - 1.0)
            lm = x1 * l1 + x2 * l2
            ln1 = (np.log(phi1 / x1) + 5.0 * q1 * np.log(th1 / phi1) + l1 - phi1 / x1 * lm
                   - q1 * np.log(th1 + th2 * t21) + q1
                   - q1 * (th1 / (th1 + th2 * t21) + th2 * t12 / (th2 + th1 * t12)))
            ln2 = (np.log(phi2 / x2) + 5.0 * q2 * np.log(th2 / phi2) + l2 - phi2 / x2 * lm
                   - q2 * np.log(th2 + th1 * t12) + q2
                   - q2 * (th2 / (th2 + th1 * t12) + th1 * t21 / (th1 + th2 * t21)))
            return np.exp(ln1), np.exp(ln2)
        raise ValueError(f"Modèle inconnu : {self.modele}")

    def K(self, x, T):
        g1, g2 = self.gammas(x, T)
        p1, p2 = self.psat(T)
        return g1 * p1 / self.P, g2 * p2 / self.P

    # ------------------------------------------------------------------ équilibre
    def point_bulle(self, x):
        """T de bulle (°C) et y pour un liquide de composition x."""
        def f(T):
            k1, k2 = self.K(x, T)
            return x * k1 + (1 - x) * k2 - 1.0
        T = brentq(f, 0.0, 200.0, xtol=1e-10)
        return T, float(x * self.K(x, T)[0])

    def _construire_tables(self):
        # azéotrope éventuel : changement de signe de y - x
        xg = np.linspace(0.02, 0.98, 49)
        ec = np.array([self.point_bulle(x)[1] - x for x in xg])
        self.x_az = self.T_az = None
        idx = np.where(np.sign(ec[:-1]) * np.sign(ec[1:]) < 0)[0]
        if len(idx):
            i = idx[0]
            self.x_az = brentq(lambda x: self.point_bulle(x)[1] - x, xg[i], xg[i + 1], xtol=1e-12)
            self.T_az = self.point_bulle(self.x_az)[0]
        xmax = self.x_az if self.x_az else 1.0
        self.x_max = xmax
        xs = np.concatenate([np.linspace(0, 0.1, 80), np.linspace(0.1, xmax, 300)[1:]])
        res = [self.point_bulle(x) if 0 < x < 1 else None for x in xs]
        T = np.array([self.cB["tnb"] - 273.15 if r is None and x == 0 else
                      (self.cA["tnb"] - 273.15 if r is None else r[0]) for x, r in zip(xs, res)])
        y = np.array([0.0 if x == 0 else (1.0 if r is None else r[1]) for x, r in zip(xs, res)])
        # températures d'ébullition exactes des corps purs (aux bornes)
        T[0] = brentq(lambda t: self.psat(t)[1] - self.P, 0, 200)
        if not self.x_az:
            T[-1] = brentq(lambda t: self.psat(t)[0] - self.P, 0, 200); y[-1] = 1.0
        else:
            y[-1] = xmax
        self._xs, self._ys, self._Ts = xs, y, T
        self._fy = PchipInterpolator(xs, y)
        self._fx = PchipInterpolator(y, xs)
        self._fT = PchipInterpolator(xs, T)

    def y_eq(self, x):    return float(self._fy(x))
    def x_eq(self, y):    return float(self._fx(y))
    def T_bulle(self, x): return float(self._fT(x))
    def T_rosee(self, y): return self.T_bulle(self.x_eq(y))

    # ------------------------------------------------------------------ enthalpies
    def _lambda(self, T, c):
        TK = np.asarray(T, dtype=float) + 273.15
        tr = np.clip((c["tc"] - TK) / (c["tc"] - c["tnb"]), 1e-6, None)
        return c["lv_nb"] * tr**0.38            # corrélation de Watson

    def hE(self, x, T):
        """Enthalpie molaire d'excès (kJ/mol) : hE = -R T² d(gE/RT)/dT, dérivée numérique."""
        def gE_RT(t):
            g1, g2 = self.gammas(x, t)
            xx = np.clip(np.asarray(x, dtype=float), 1e-12, 1 - 1e-12)
            return xx * np.log(g1) + (1 - xx) * np.log(g2)
        dT = 0.05
        d = (gE_RT(np.asarray(T) + dT) - gE_RT(np.asarray(T) - dT)) / (2 * dT)
        TK = np.asarray(T, dtype=float) + 273.15
        return -R_GAZ * TK**2 * d

    def h_liq(self, x, T, avec_hE=True):
        x = np.asarray(x, dtype=float)
        h = (x * self.cA["cp_l"] + (1 - x) * self.cB["cp_l"]) * (np.asarray(T) - T_REF)
        return h + self.hE(x, T) if avec_hE else h

    def H_vap(self, y, T):
        y = np.asarray(y, dtype=float)
        T = np.asarray(T, dtype=float)
        hA = self.cA["cp_l"] * (T - T_REF) + self._lambda(T, self.cA)
        hB = self.cB["cp_l"] * (T - T_REF) + self._lambda(T, self.cB)
        return y * hA + (1 - y) * hB

    # enthalpies sur les courbes de saturation (pour Ponchon-Savarit)
    def h_sat(self, x): return float(self.h_liq(x, self.T_bulle(x)))
    def H_sat(self, y): return float(self.H_vap(y, self.T_rosee(y)))

    def enthalpie_alimentation(self, zF, q):
        """q = 1 : liquide saturé ; q = 0 : vapeur saturée ; 0<q<1 : mélange."""
        hl = self.h_sat(zF)
        Hv = float(self.H_vap(zF, self.T_rosee(zF)))
        return hl + (1.0 - q) * (Hv - hl)
