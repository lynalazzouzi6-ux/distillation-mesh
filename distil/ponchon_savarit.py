"""
Méthode de Ponchon-Savarit (binaire, diagramme enthalpie-composition), sur le même modèle
thermodynamique que le MESH (NRTL/UNIQUAC + enthalpies avec hE) : les deux méthodes doivent
donc donner la même solution (test d'équivalence dans validation.py).
"""
from scipy.optimize import brentq


def construire(sys, F, zF, xD, xW, R, q, n_max=100, n_alim_force=None):
    """Construit les plateaux depuis le distillat. Retourne (resultats, plateaux, droites_op)."""
    xm = sys.x_max
    if xD >= xm:
        raise ValueError(f"xD={xD} >= limite de l'azéotrope ({xm:.3f})")
    D = F * (zF - xW) / (xD - xW)
    W = F - D
    hF = sys.enthalpie_alimentation(zF, q)
    hD, H1 = sys.h_sat(xD), sys.H_sat(xD)
    H_DELTA = hD + (R + 1) * (H1 - hD)
    H_DELTA_P = H_DELTA + (hF - H_DELTA) * (xW - xD) / (zF - xD)
    ligne_DF = lambda x: H_DELTA + (hF - H_DELTA) * (x - xD) / (zF - xD)

    pole, apres, y = (xD, H_DELTA), False, xD
    plateaux, ops, n_alim = [], [], None
    for n in range(1, n_max + 1):
        x = sys.x_eq(y)
        h, H = sys.h_sat(x), sys.H_sat(y)
        if n_alim_force is None:
            bascule = (h - ligne_DF(x)) * (H - ligne_DF(y)) <= 0
        else:
            bascule = (n == n_alim_force)
        if (not apres) and bascule:
            apres, pole, n_alim = True, (xW, H_DELTA_P), n
        plateaux.append((x, h, y, H))
        if x <= xW + 1e-7:
            break
        px, pH = pole
        pente = (h - pH) / (x - px)
        f = lambda yy: sys.H_sat(yy) - (pH + pente * (yy - px))
        lo, hi = x + 1e-9, xm - 1e-9
        if f(lo) * f(hi) >= 0:
            raise RuntimeError(f"Pincement (R={R:.3f}, étage {n})")
        y_new = brentq(f, lo, hi)
        ops.append(((x, h), (y_new, sys.H_sat(y_new)), pole))
        y = y_new
    else:
        raise RuntimeError("n_max atteint : augmente R")

    nt = len(plateaux)
    frac = (plateaux[-2][0] - xW) / (plateaux[-2][0] - plateaux[-1][0]) if nt >= 2 else 1.0
    hW = sys.h_sat(xW)
    res = dict(D=D, W=W, hF=hF, H_DELTA=H_DELTA, H_DELTA_P=H_DELTA_P,
               Qc=D * (H_DELTA - hD), Qr=W * (hW - H_DELTA_P),          # MJ/h
               n_etages=nt, n_frac=nt - 1 + frac, n_alim=n_alim,
               x_dernier=plateaux[-1][0], xD=xD, xW=xW, zF=zF, R=R, q=q)
    return res, plateaux, ops


def r_min(sys, F, zF, xD, xW, q, tol=1e-3):
    """Reflux minimal (dichotomie sur la faisabilité avec un grand nombre d'étages)."""
    def ok(R):
        try:
            construire(sys, F, zF, xD, xW, R, q, n_max=400); return True
        except Exception:
            return False
    lo, hi = 0.02, 30.0
    if not ok(hi):
        raise RuntimeError("Aucun reflux fini ne convient (xD trop proche de l'azéotrope ?)")
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if ok(mid): hi = mid
        else: lo = mid
    return hi
