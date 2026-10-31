"""
Validation du modèle :
  1. VLE : comparaison à des données expérimentales (si disponibles) et à la bibliothèque `thermo` (si installée)
  2. Équivalence Ponchon-Savarit <-> MESH (même modèle thermodynamique : les profils doivent coïncider)
  3. Fermeture des bilans (matière, composant, énergie) du MESH
"""
import numpy as np
from . import ponchon_savarit as ps, mesh


def erreur_vle(sys):
    """RMSE sur y et T par rapport aux données expérimentales (None si pas de données)."""
    if not sys.donnees_exp:
        return None
    d = np.array(sys.donnees_exp)
    pts = [sys.point_bulle(x) for x in d[:, 0]]
    T = np.array([p[0] for p in pts]); y = np.array([p[1] for p in pts])
    return dict(rmse_y=float(np.sqrt(np.mean((y - d[:, 1]) ** 2))),
                rmse_T=float(np.sqrt(np.mean((T - d[:, 2]) ** 2))),
                x=d[:, 0], y_exp=d[:, 1], T_exp=d[:, 2], y_calc=y, T_calc=T,
                n=len(d), x_az=sys.x_az, T_az=sys.T_az)


def comparaison_thermo(sys, x=0.3, T=70.0):
    """Compare nos γ et hE avec la bibliothèque thermo (optionnelle). Retourne None si absente."""
    try:
        from thermo.nrtl import NRTL
        from thermo.uniquac import UNIQUAC
    except Exception:
        return None
    try:
        if sys.modele == "NRTL":
            p = sys.params
            m = NRTL(T=T + 273.15, xs=[x, 1 - x], tau_bs=[[0, p["b12"]], [p["b21"], 0]],
                     alpha_cs=[[0, p["alpha"]], [p["alpha"], 0]])
        else:
            p = sys.params
            m = UNIQUAC(T=T + 273.15, xs=[x, 1 - x], rs=[sys.cA["r"], sys.cB["r"]],
                        qs=[sys.cA["q"], sys.cB["q"]], tau_bs=[[0, p["b12"]], [p["b21"], 0]])
        g_ref = m.gammas(); g = sys.gammas(x, T)
        return dict(gamma_ref=[float(v) for v in g_ref], gamma_nous=[float(g[0]), float(g[1])],
                    ecart_max=float(max(abs(g_ref[0] - g[0]), abs(g_ref[1] - g[1]))),
                    hE_ref=float(m.HE() / 1000), hE_nous=float(sys.hE(x, T)))
    except Exception:
        return None


def equivalence_ps_mesh(sys, F, zF, xD, xW, q, R, N, NF):
    """
    Résout le MESH pour (N, NF, R), puis refait Ponchon-Savarit avec les puretés obtenues.
    PS doit trouver exactement N étages, l'alimentation à NF et les mêmes compositions.
    """
    D = F * (zF - xW) / (xD - xW)
    sol = mesh.resoudre(sys, N, NF, F, zF, q, D, R, xD_guess=xD, xW_guess=xW)
    if not sol.convergee:
        return None
    res, pl, ops = ps.construire(sys, F, zF, sol.xD, sol.xW, R, q, n_alim_force=NF)
    xps = np.array([p[0] for p in pl])
    n = min(len(xps), N)
    return dict(sol=sol, res_ps=res, plateaux=pl, n_frac_ps=res["n_frac"], N=N,
                ecart_x=float(np.max(np.abs(sol.x[:n] - xps[:n]))),
                ecart_Qc=abs(sol.Qc - res["Qc"]) / sol.Qc, ecart_Qr=abs(sol.Qr - res["Qr"]) / sol.Qr)
