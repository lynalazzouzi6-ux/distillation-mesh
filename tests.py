"""
Tests automatiques : python tests.py
Vérifie que le projet est correct (thermodynamique, MESH, équivalence avec Ponchon-Savarit, optimisation).
"""
import sys as _sys
try: _sys.stdout.reconfigure(encoding="utf-8")
except Exception: pass
import numpy as np
from distil.systemes import charger, SYSTEMES
from distil import ponchon_savarit as ps, mesh, optimisation as opt, validation

ok = 0
def verifie(nom, condition, detail=""):
    global ok
    print(("  [OK]  " if condition else "  [ÉCHEC] ") + nom + (f"  ({detail})" if detail else ""))
    if not condition: raise SystemExit(1)
    ok += 1

print("1. Thermodynamique")
s = charger("ethanol_eau")
verifie("azéotrope éthanol/eau entre 0,85 et 0,90 (exp. : 0,894)", 0.85 < s.x_az < 0.90, f"x_az = {s.x_az:.4f}")
verifie("T azéotrope proche de 78,15 °C (exp.)", abs(s.T_az - 78.15) < 0.3, f"{s.T_az:.2f} °C")
verifie("T de bulle du pur eau ≈ 100 °C", abs(s.T_bulle(0.0) - 100.0) < 0.2, f"{s.T_bulle(0.0):.2f} °C")
v = validation.erreur_vle(s)
verifie("RMSE y vs données exp. < 0,02", v["rmse_y"] < 0.02, f"{v['rmse_y']:.4f}")
verifie("RMSE T vs données exp. < 0,5 °C", v["rmse_T"] < 0.5, f"{v['rmse_T']:.2f}")
c = validation.comparaison_thermo(s)
if c: verifie("γ identiques à la bibliothèque thermo (NRTL)", c["ecart_max"] < 1e-8, f"{c['ecart_max']:.1e}")
else: print("  (thermo non installé : comparaison ignorée)")
su = charger("ethanol_eau_uniquac")
c = validation.comparaison_thermo(su)
if c: verifie("γ identiques à la bibliothèque thermo (UNIQUAC)", c["ecart_max"] < 1e-8, f"{c['ecart_max']:.1e}")
x = 0.4; T = 80.0
verifie("Gibbs-Duhem : x1 dlnγ1 + x2 dlnγ2 = 0",
        abs(x * (np.log(s.gammas(x + 1e-6, T)[0]) - np.log(s.gammas(x - 1e-6, T)[0])) / 2e-6
            + (1 - x) * (np.log(s.gammas(x + 1e-6, T)[1]) - np.log(s.gammas(x - 1e-6, T)[1])) / 2e-6) < 1e-6)

print("2. MESH et bilans")
d = s.defauts; F, zF, xD, xW, q = d["F"], d["zF"], d["xD"], d["xW"], d["q"]
D = F * (zF - xW) / (xD - xW)
sol = mesh.resoudre(s, 19, 17, F, zF, q, D, 1.3, xD_guess=xD, xW_guess=xW)
verifie("MESH converge", sol.convergee, f"{sol.iterations} itérations")
verifie("bilan matière fermé", sol.bilans["matiere"] < 1e-9)
verifie("bilan composant fermé", sol.bilans["composant"] < 1e-9)
verifie("bilan énergie fermé", sol.bilans["energie"] < 1e-9)
verifie("Qc et Qr positifs et du même ordre", sol.Qc > 0 and sol.Qr > 0 and abs(sol.Qr / sol.Qc - 1) < 0.2)
verifie("profil x décroissant de la tête au pied", np.all(np.diff(sol.x) < 0))
verifie("profil T croissant de la tête au pied", np.all(np.diff(sol.T) > 0))

print("3. Équivalence MESH <-> Ponchon-Savarit (tous les systèmes)")
for k in SYSTEMES:
    sy = charger(k); dd = sy.defauts
    rm = ps.r_min(sy, dd["F"], dd["zF"], dd["xD"], dd["xW"], dd["q"])
    r0, _, _ = ps.construire(sy, dd["F"], dd["zF"], dd["xD"], dd["xW"], 1.5 * rm, dd["q"])
    eq = validation.equivalence_ps_mesh(sy, dd["F"], dd["zF"], dd["xD"], dd["xW"], dd["q"], 1.5 * rm, r0["n_etages"], r0["n_alim"])
    verifie(f"{k}: même solution (écart x < 1e-5)", eq is not None and eq["ecart_x"] < 1e-5 and abs(eq["n_frac_ps"] - eq["N"]) < 1e-3,
            f"écart x = {eq['ecart_x']:.1e}")

print("4. Optimisation")
rm = ps.r_min(s, F, zF, xD, xW, q)
r0, _, _ = ps.construire(s, F, zF, xD, xW, 1.5 * rm, q)
C = opt.Concepteur(s, F, zF, xD, xW, q)
designs = C.balayage(r0["n_frac"], rm, N_max_extra=22, patience=8)
best = min(designs, key=lambda dd: dd["cout"]["tac"])
verifie("l'optimum n'est pas au bord du balayage", designs[0]["cout"]["N"] < best["cout"]["N"] < designs[-1]["cout"]["N"], f"N* = {best['cout']['N']}")
verifie("R optimal entre 1,05 et 1,6 × R_min", 1.05 < best["cout"]["R"] / rm < 1.6, f"{best['cout']['R']/rm:.2f}")
verifie("pureté respectée", best["sol"].xD >= xD - 1e-6 and best["sol"].xW <= xW + 1e-4)
v = C.nlp_verification(best["sol"])
verifie("le NLP (SLSQP) confirme l'optimum", abs(v["tac"] - best["cout"]["tac"]) / best["cout"]["tac"] < 1e-3)
print(f"\nTous les tests sont passés ({ok}).")
