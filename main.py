"""
PROJET : simulateur rigoureux de distillation binaire non idéale (MESH) avec optimisation CAPEX/OPEX.

Utilisation (depuis le dossier du projet) :
    python main.py                          # éthanol/eau (par défaut)
    python main.py --systeme methanol_eau
    python main.py --systeme tous           # les 4 systèmes
    python main.py --rapide                 # sans l'étude de sensibilité (plus court)
    python main.py --xD 0.85 --zF 0.25 --q 0.5 --vapeur 40

Résultats : dossier  resultats/<systeme>/  (figures PNG, rapport.md, designs.csv)
"""
import argparse, csv, os, sys as _sys, time

for _flux in (_sys.stdout, _sys.stderr):          # affichage correct des accents/symboles sous Windows
    try: _flux.reconfigure(encoding="utf-8")
    except Exception: pass
from dataclasses import replace
import numpy as np

from distil.systemes import charger, SYSTEMES
from distil import ponchon_savarit as ps, optimisation as opt, economics, validation, graphiques as gr


def executer(cle, args):
    t0 = time.time()
    sys = charger(cle)
    d = dict(sys.defauts)
    for k_arg, k in (("F", "F"), ("zF", "zF"), ("xD", "xD"), ("xW", "xW"), ("q", "q")):
        if getattr(args, k_arg) is not None: d[k] = getattr(args, k_arg)
    F, zF, xD, xW, q = d["F"], d["zF"], d["xD"], d["xW"], d["q"]
    eco = economics.ParamsEco(prix_vapeur=args.vapeur)
    dossier = os.path.join(args.sortie, cle)
    os.makedirs(dossier, exist_ok=True)
    chemin = lambda nom: os.path.join(dossier, nom)
    rap = []
    def pr(s=""):
        print(s); rap.append(s)

    pr(f"# Distillation {sys.nomA}/{sys.nomB} — {sys.description}")
    pr(f"Modèle : {sys.modele} | P = {sys.P/760:.2f} atm | F = {F} kmol/h, zF = {zF}, q = {q}")
    pr(f"Spécifications : xD ≥ {xD}, xW ≤ {xW}")
    if sys.x_az:
        pr(f"Azéotrope du modèle : x = {sys.x_az:.4f}, T = {sys.T_az:.2f} °C  (xD doit rester en dessous)")
    else:
        pr("Pas d'azéotrope : séparation possible jusqu'à la pureté complète.")
    if sys.x_az and xD >= sys.x_az:
        pr("ERREUR : xD dépasse l'azéotrope, change la spécification."); return None

    # ------------------------------------------------------------------ 1. VLE
    pr("\n## 1. Validation de l'équilibre liquide-vapeur")
    val = validation.erreur_vle(sys)
    if val:
        pr(f"Écart avec {val['n']} points expérimentaux (données de manuel) : RMSE(y) = {val['rmse_y']:.4f}, RMSE(T) = {val['rmse_T']:.2f} °C")
    else:
        pr("Pas de données expérimentales enregistrées pour ce système (ajoute-les dans systemes.py).")
    cmp_t = validation.comparaison_thermo(sys)
    if cmp_t:
        pr(f"Comparaison avec la bibliothèque thermo (x=0.3, 70 °C) : γ nous = {cmp_t['gamma_nous'][0]:.5f}/{cmp_t['gamma_nous'][1]:.5f}, "
           f"thermo = {cmp_t['gamma_ref'][0]:.5f}/{cmp_t['gamma_ref'][1]:.5f} (écart max {cmp_t['ecart_max']:.1e}) ; "
           f"hE nous = {cmp_t['hE_nous']:.4f}, thermo = {cmp_t['hE_ref']:.4f} kJ/mol")
    gr.fig_vle(sys, chemin("01_equilibre_VLE.png"), val)

    # ------------------------------------------------------------------ 2. Ponchon-Savarit
    pr("\n## 2. Ponchon-Savarit")
    Rm = ps.r_min(sys, F, zF, xD, xW, q)
    R_base = 1.5 * Rm
    res, pl, ops = ps.construire(sys, F, zF, xD, xW, R_base, q)
    pr(f"R_min = {Rm:.3f} | à R = 1,5 R_min = {R_base:.3f} : {res['n_frac']:.2f} étages théoriques, alimentation étage {res['n_alim']}, "
       f"Qc = {res['Qc']/3.6:.0f} kW, Qr = {res['Qr']/3.6:.0f} kW")
    gr.fig_ponchon(sys, res, pl, ops, f"Ponchon-Savarit : R = {R_base:.2f} ({res['n_frac']:.1f} étages)", chemin("02_ponchon_savarit.png"))

    courbes = []
    for qq in (1.0, 0.5, 0.0):
        try:
            rmq = ps.r_min(sys, F, zF, xD, xW, qq)
            Rs = np.linspace(1.03 * rmq, max(3.5, 2.5 * rmq), 30)
            Ns, Qrs = [], []
            for r in Rs:
                try:
                    rr, _, _ = ps.construire(sys, F, zF, xD, xW, r, qq, n_max=150)
                    Ns.append(rr["n_frac"]); Qrs.append(rr["Qr"] / 3.6)
                except Exception:
                    Ns.append(np.nan); Qrs.append(np.nan)
            courbes.append((qq, rmq, Rs, Ns, Qrs))
            pr(f"  q = {qq}: R_min = {rmq:.3f}")
        except Exception as e:
            pr(f"  q = {qq}: non calculable ({e})")
    gr.fig_effet_R_q(courbes, chemin("03_effet_R_et_q.png"))

    # ------------------------------------------------------------------ 3. MESH vs PS
    pr("\n## 3. MESH rigoureux et équivalence avec Ponchon-Savarit")
    eq = validation.equivalence_ps_mesh(sys, F, zF, xD, xW, q, R_base, res["n_etages"], res["n_alim"])
    if eq:
        s = eq["sol"]
        pr(f"MESH à N = {s.N}, N_F = {s.NF}, R = {R_base:.3f} : xD = {s.xD:.5f}, xW = {s.xW:.5f}, convergé en {s.iterations} itérations "
           f"(résidu {s.residu:.1e})")
        pr(f"Ponchon-Savarit refait avec ces puretés : {eq['n_frac_ps']:.3f} étages (attendu {eq['N']}) ; écart max sur x = {eq['ecart_x']:.1e}, "
           f"écart relatif Qc = {eq['ecart_Qc']:.1e}, Qr = {eq['ecart_Qr']:.1e}")
        pr(f"Fermeture des bilans : matière {s.bilans['matiere']:.1e}, composant {s.bilans['composant']:.1e}, énergie (relative) {s.bilans['energie']:.1e}")
        gr.fig_equivalence(eq, sys, chemin("04_equivalence_PS_MESH.png"))
    else:
        pr("MESH non convergé pour ce cas de base.")

    # ------------------------------------------------------------------ 4. optimisation
    pr("\n## 4. Optimisation CAPEX + OPEX (MESH, variables N, N_F, R)")
    C = opt.Concepteur(sys, F, zF, xD, xW, q, eco)
    designs = C.balayage(res["n_frac"], Rm, N_max_extra=(14 if args.rapide else 30), patience=(7 if args.rapide else 100))
    tacs = [dd["cout"]["tac"] for dd in designs]
    ib = int(np.argmin(tacs)); best = designs[ib]; cb = best["cout"]; sb = best["sol"]
    pr(f"{len(designs)} conceptions évaluées ({C.n_solves} résolutions MESH).")
    pr(f"**Optimum : N = {sb.N} étages théoriques ({cb['n_reels']} plateaux réels, E_o = {eco.E_o}), alimentation étage {sb.NF}, "
       f"R = {sb.R:.3f} ({sb.R/Rm:.2f} × R_min)**")
    pr(f"Colonne : Ø {cb['diam']:.2f} m, H {cb['haut']:.1f} m | Qc = {cb['Qc']:.0f} kW, Qr = {cb['Qr']:.0f} kW | xD = {sb.xD:.4f}, xW = {sb.xW:.4f}")
    pr(f"CAPEX = {cb['capex']/1e3:.0f} k$ ({cb['capex_annuel']/1e3:.0f} k$/an) | OPEX = {cb['opex']/1e3:.0f} k$/an | TAC = {cb['tac']/1e3:.0f} k$/an")
    res_o, _, _ = ps.construire(sys, F, zF, xD, xW, sb.R, q)
    pr(f"Comparaison au même R : Ponchon-Savarit donne {res_o['n_frac']:.2f} étages (MESH : {sb.N} entiers).")
    try:
        v = C.nlp_verification(sb)
        pr(f"Vérification NLP (SLSQP sur R et D, contraintes xD ≥ {xD}, xW ≤ {xW}) : R = {v['R']:.4f}, D = {v['D']:.3f} kmol/h, "
           f"TAC = {v['tac']/1e3:.1f} k$/an, succès = {v['succes']} -> {'optimum confirmé' if abs(v['tac']-cb['tac'])/cb['tac'] < 1e-3 else 'écart à examiner'}")
    except Exception as e:
        pr(f"Vérification NLP non exécutée ({e})")
    gr.fig_optimisation(designs, ib, eco, chemin("05_optimisation.png"))
    gr.fig_profils(sb, sys, f"Solution optimale : N = {sb.N}, N_F = {sb.NF}, R = {sb.R:.2f}", chemin("06_profils_colonne.png"))
    rr, pl2, ops2 = ps.construire(sys, F, zF, xD, xW, sb.R, q)
    gr.fig_ponchon(sys, rr, pl2, ops2, f"Solution optimale (R = {sb.R:.2f}) sur le diagramme H-x-y", chemin("07_ponchon_solution.png"))

    with open(chemin("designs.csv"), "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f, delimiter=";")
        w.writerow(["N", "NF", "R", "plateaux_reels", "diametre_m", "hauteur_m", "Qc_kW", "Qr_kW", "CAPEX_k$", "OPEX_k$/an", "TAC_k$/an"])
        for dd in designs:
            c = dd["cout"]
            w.writerow([c["N"], c["NF"], f"{c['R']:.4f}", c["n_reels"], f"{c['diam']:.3f}", f"{c['haut']:.2f}",
                        f"{c['Qc']:.1f}", f"{c['Qr']:.1f}", f"{c['capex']/1e3:.1f}", f"{c['opex']/1e3:.1f}", f"{c['tac']/1e3:.1f}"])

    # ------------------------------------------------------------------ 5. sensibilité
    if not args.rapide:
        pr("\n## 5. Sensibilité au prix de la vapeur (mêmes conceptions, coûts recalculés)")
        prix_l, nopt_l, rr_l, tac_l = [], [], [], []
        for p in (5, 10, 25, 50, 100, 200):
            e2 = replace(eco, prix_vapeur=p)
            cs = [economics.evaluer(dd["sol"], e2) for dd in designs]
            k = int(np.argmin([c["tac"] for c in cs]))
            prix_l.append(p); nopt_l.append(cs[k]["N"]); rr_l.append(cs[k]["R"] / Rm); tac_l.append(cs[k]["tac"])
            bord = "  (bord du balayage : l'optimum réel est peut-être plus loin)" if cs[k]["N"] >= max(c["N"] for c in cs) else ""
            pr(f"  vapeur {p:4d} $/t -> N_opt = {cs[k]['N']}, N_F = {cs[k]['NF']}, R_opt = {cs[k]['R']:.2f} ({cs[k]['R']/Rm:.2f} × R_min), TAC = {cs[k]['tac']/1e3:.0f} k$/an{bord}")
        gr.fig_sensibilite(prix_l, rr_l, nopt_l, tac_l, chemin("08_sensibilite_vapeur.png"))

    pr(f"\nTemps de calcul : {time.time()-t0:.0f} s")
    pr("\n## Limites du modèle")
    pr("- Étages théoriques, P constante (pas de perte de charge), condenseur total, une alimentation, mélange binaire.")
    pr("- Cp liquides constants, chaleurs latentes par la corrélation de Watson, vapeur idéale, enthalpie d'excès via le modèle de γ.")
    pr("- Paramètres NRTL/UNIQUAC de type ChemSep (à remplacer par ceux de ton cours si besoin) ; coûts = ordres de grandeur.")
    with open(chemin("rapport.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(rap))
    return dict(cle=cle, N=sb.N, NF=sb.NF, R=sb.R, Rm=Rm, tac=cb["tac"])


def main():
    ap = argparse.ArgumentParser(description="Simulateur MESH de distillation binaire + optimisation")
    ap.add_argument("--systeme", default="ethanol_eau", help=f"{', '.join(SYSTEMES)} ou 'tous'")
    ap.add_argument("--F", type=float); ap.add_argument("--zF", type=float)
    ap.add_argument("--xD", type=float); ap.add_argument("--xW", type=float); ap.add_argument("--q", type=float)
    ap.add_argument("--vapeur", type=float, default=25.0, help="prix de la vapeur ($/t)")
    ap.add_argument("--rapide", action="store_true", help="version courte (sans sensibilité)")
    ap.add_argument("--sortie", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "resultats"))
    args = ap.parse_args()
    cles = list(SYSTEMES) if args.systeme == "tous" else [args.systeme]
    synth = []
    for k in cles:
        print("=" * 78)
        r = executer(k, args)
        if r: synth.append(r)
        print()
    if len(synth) > 1:
        print("=" * 78 + "\nSYNTHÈSE")
        for r in synth:
            print(f"  {r['cle']:22s} N = {r['N']:2d}, N_F = {r['NF']:2d}, R = {r['R']:.2f} ({r['R']/r['Rm']:.2f} × R_min), TAC = {r['tac']/1e3:.0f} k$/an")
    print(f"\nTerminé. Résultats dans : {args.sortie}")


if __name__ == "__main__":
    main()
