"""
Optimisation de la conception sur le modèle MESH rigoureux.
Variables de décision : N (étages théoriques), N_F (étage d'alimentation), R (taux de reflux).
Contraintes de pureté : xD >= xD_spec et xW <= xW_spec.
Méthode :
  - D est fixé par les bilans globaux avec les deux puretés spécifiées ;
  - pour (N, N_F) donnés, le plus petit R qui donne xD = xD_spec est trouvé par recherche de racine
    (Brent) sur le MESH complet ;
  - N_F est optimisé par descente discrète (celui qui rend xD maximal à R donné) ;
  - N est parcouru (optimisation discrète) et le coût annuel total TAC = OPEX + CAPEX annualisé
    est minimisé ;
  - un NLP continu (SLSQP, variables R et D, contraintes de pureté) vérifie l'optimalité locale.
"""
import numpy as np
from scipy.optimize import brentq, minimize
from . import mesh, economics


class Concepteur:
    def __init__(self, sys, F, zF, xD_spec, xW_spec, q, eco=None):
        self.s, self.F, self.zF, self.xDs, self.xWs, self.q = sys, F, zF, xD_spec, xW_spec, q
        self.D = F * (zF - xW_spec) / (xD_spec - xW_spec)
        self.eco = eco or economics.ParamsEco()
        self._cache = {}           # (N, NF) -> (R, u)
        self.n_solves = 0

    # ------------------------------------------------------------------ un calcul MESH
    def solve(self, N, NF, R, D=None):
        D = self.D if D is None else D
        key = (N, NF)
        u0 = None
        if key in self._cache:
            R_old, u_old = self._cache[key]
            u0 = u_old.copy()
            dR = (R - R_old) * self.D
            u0[2 * N:3 * N - 1] += dR              # L_1..L_{N-1}
            u0[3 * N:4 * N - 1] += dR              # V_2..V_N
        self.n_solves += 1
        sol = mesh.resoudre(self.s, N, NF, self.F, self.zF, self.q, D, R, u0=u0,
                            xD_guess=self.xDs, xW_guess=self.xWs)
        if not sol.convergee and u0 is not None:   # repli : estimation initiale
            sol = mesh.resoudre(self.s, N, NF, self.F, self.zF, self.q, D, R,
                                xD_guess=self.xDs, xW_guess=self.xWs)
        if sol.convergee and D == self.D:
            self._cache[key] = (R, sol.u)
        return sol

    def xD(self, N, NF, R):
        sol = self.solve(N, NF, R)
        return sol.xD if sol.convergee else np.nan

    # ------------------------------------------------------------------ R requis pour (N, NF)
    def R_requis(self, N, NF, R_lo, R_hi):
        f = lambda R: self.xD(N, NF, R) - self.xDs
        # borne haute : on l'élargit progressivement (évite les R extrêmes où le MESH diverge)
        hi = None
        for cand in sorted({min(R_hi, m * R_lo / 0.3) for m in (1.0, 1.6, 2.5, 4.0, 8.0)} | {R_hi}):
            v = f(cand)
            if np.isfinite(v) and v >= 0:
                hi = cand; break
        if hi is None:
            return None                              # N trop petit : pureté inatteignable
        lo = min(R_lo, 0.9 * hi)
        flo = f(lo)
        while np.isfinite(flo) and flo > 0 and lo > 0.02:
            lo *= 0.5; flo = f(lo)
        if not np.isfinite(flo):
            return None
        g = lambda R: (lambda v: v if np.isfinite(v) else 1.0)(f(R))
        R = brentq(g, lo, hi, xtol=1e-6, rtol=1e-8)
        sol = self.solve(N, NF, R)
        return sol if sol.convergee else None

    # ------------------------------------------------------------------ N_F optimal pour (N, R)
    def meilleur_NF(self, N, R, NF0):
        NF = int(np.clip(NF0, 2, N - 1))
        best = self.xD(N, NF, R)
        for sens in (+1, -1):
            while True:
                cand = NF + sens
                if not (2 <= cand <= N - 1): break
                v = self.xD(N, cand, R)
                if np.isfinite(v) and (not np.isfinite(best) or v > best + 1e-9):
                    NF, best = cand, v
                else:
                    break
        return NF

    # ------------------------------------------------------------------ une conception pour N
    def concevoir(self, N, R_ref, NF0, R_hi):
        NF = self.meilleur_NF(N, R_ref, NF0)
        sol = self.R_requis(N, NF, 0.3 * R_ref, R_hi)
        if sol is None:
            return None
        NF2 = self.meilleur_NF(N, sol.R, NF)         # ré-optimisation de l'alimentation à R requis
        if NF2 != NF:
            s2 = self.R_requis(N, NF2, 0.5 * sol.R, 1.5 * sol.R + 0.5)
            if s2 is not None and s2.R < sol.R:
                sol = s2
        return sol

    # ------------------------------------------------------------------ balayage en N
    def balayage(self, N_ps, R_min_ps, N_max_extra=28, patience=7, verbose=False):
        R_hi = 25.0 * max(R_min_ps, 0.3)
        designs = []
        # plus petit N faisable
        N = max(4, int(0.45 * N_ps))
        sol = None
        while N < 4 * N_ps + 20:
            NF0 = max(2, int(round(0.7 * N)))
            sol = self.concevoir(N, 1.5 * R_min_ps, NF0, R_hi)
            if sol is not None: break
            N += 1
        if sol is None:
            raise RuntimeError("Aucun N faisable trouvé")
        N_first = N
        best_tac, sans_amelioration = np.inf, 0
        R_ref, NF_ref = sol.R, sol.NF
        while N <= N_first + N_max_extra:
            if sol is None:
                sol = self.concevoir(N, R_ref, NF_ref, R_hi)
            if sol is None:
                N += 1; continue
            c = economics.evaluer(sol, self.eco)
            designs.append(dict(sol=sol, cout=c))
            if verbose:
                print(f"  N={N:3d} NF={sol.NF:3d} R={sol.R:.3f} TAC={c['tac']/1e3:.1f} k$/an")
            if c["tac"] < best_tac - 1e-6:
                best_tac, sans_amelioration = c["tac"], 0
            else:
                sans_amelioration += 1
                if sans_amelioration >= patience: break
            R_ref, NF_ref = sol.R, sol.NF
            N += 1; sol = None
        return designs

    # ------------------------------------------------------------------ NLP de vérification
    def nlp_verification(self, sol0):
        """SLSQP sur (R, D) à (N, NF) fixés : min TAC s.t. xD >= xD_spec, xW <= xW_spec."""
        N, NF = sol0.N, sol0.NF
        def evalue(v):
            R, D = v
            s = self.solve(N, NF, R, D=D)
            return s if s.convergee else None
        def obj(v):
            s = evalue(v)
            return economics.evaluer(s, self.eco)["tac"] / 1e5 if s else 1e6
        def c_xD(v):
            s = evalue(v); return (s.xD - self.xDs) * 100 if s else -1.0
        def c_xW(v):
            s = evalue(v); return (self.xWs - s.xW) * 100 if s else -1.0
        v0 = np.array([sol0.R, self.D])
        r = minimize(obj, v0, method="SLSQP",
                     constraints=[dict(type="ineq", fun=c_xD), dict(type="ineq", fun=c_xW)],
                     bounds=[(0.05, 40.0), (0.5 * self.D, 1.5 * self.D)],
                     options=dict(ftol=1e-9, maxiter=60, eps=1e-6))
        s = evalue(r.x)
        return dict(R=float(r.x[0]), D=float(r.x[1]), tac=obj(r.x) * 1e5, succes=bool(r.success),
                    xD=s.xD if s else np.nan, xW=s.xW if s else np.nan, message=str(r.message))
