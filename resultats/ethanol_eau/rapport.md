# Distillation Éthanol/Eau — Éthanol/eau, NRTL, azéotrope vers x = 0,88
Modèle : NRTL | P = 1.00 atm | F = 100.0 kmol/h, zF = 0.2, q = 1.0
Spécifications : xD ≥ 0.8, xW ≤ 0.01
Azéotrope du modèle : x = 0.8814, T = 78.10 °C  (xD doit rester en dessous)

## 1. Validation de l'équilibre liquide-vapeur
Écart avec 14 points expérimentaux (données de manuel) : RMSE(y) = 0.0089, RMSE(T) = 0.29 °C
Comparaison avec la bibliothèque thermo (x=0.3, 70 °C) : γ nous = 1.75984/1.20021, thermo = 1.75984/1.20021 (écart max 6.7e-16) ; hE nous = 0.6440, thermo = 0.6440 kJ/mol

## 2. Ponchon-Savarit
R_min = 0.989 | à R = 1,5 R_min = 1.483 : 14.85 étages théoriques, alimentation étage 13, Qc = 647 kW, Qr = 653 kW
  q = 1.0: R_min = 0.989
  q = 0.5: R_min = 1.695
  q = 0.0: R_min = 3.587

## 3. MESH rigoureux et équivalence avec Ponchon-Savarit
MESH à N = 15, N_F = 13, R = 1.483 : xD = 0.80108, xW = 0.00966, convergé en 5 itérations (résidu 4.9e-14)
Ponchon-Savarit refait avec ces puretés : 15.000 étages (attendu 15) ; écart max sur x = 3.8e-08, écart relatif Qc = 3.2e-10, Qr = 9.5e-10
Fermeture des bilans : matière 0.0e+00, composant 1.7e-14, énergie (relative) 3.2e-14

## 4. Optimisation CAPEX + OPEX (MESH, variables N, N_F, R)
31 conceptions évaluées (529 résolutions MESH).
**Optimum : N = 19 étages théoriques (26 plateaux réels, E_o = 0.7), alimentation étage 17, R = 1.267 (1.28 × R_min)**
Colonne : Ø 0.84 m, H 18.6 m | Qc = 590 kW, Qr = 597 kW | xD = 0.8000, xW = 0.0100
CAPEX = 492 k$ (80 k$/an) | OPEX = 222 k$/an | TAC = 302 k$/an
Comparaison au même R : Ponchon-Savarit donne 19.00 étages (MESH : 19 entiers).
Vérification NLP (SLSQP sur R et D, contraintes xD ≥ 0.8, xW ≤ 0.01) : R = 1.2673, D = 24.051 kmol/h, TAC = 302.0 k$/an, succès = True -> optimum confirmé

## 5. Sensibilité au prix de la vapeur (mêmes conceptions, coûts recalculés)
  vapeur    5 $/t -> N_opt = 15, N_F = 13, R_opt = 1.45 (1.47 × R_min), TAC = 135 k$/an
  vapeur   10 $/t -> N_opt = 15, N_F = 13, R_opt = 1.45 (1.47 × R_min), TAC = 179 k$/an
  vapeur   25 $/t -> N_opt = 19, N_F = 17, R_opt = 1.27 (1.28 × R_min), TAC = 302 k$/an
  vapeur   50 $/t -> N_opt = 24, N_F = 21, R_opt = 1.16 (1.17 × R_min), TAC = 497 k$/an
  vapeur  100 $/t -> N_opt = 29, N_F = 26, R_opt = 1.11 (1.12 × R_min), TAC = 876 k$/an
  vapeur  200 $/t -> N_opt = 36, N_F = 33, R_opt = 1.06 (1.08 × R_min), TAC = 1616 k$/an

Temps de calcul : 11 s

## Limites du modèle
- Étages théoriques, P constante (pas de perte de charge), condenseur total, une alimentation, mélange binaire.
- Cp liquides constants, chaleurs latentes par la corrélation de Watson, vapeur idéale, enthalpie d'excès via le modèle de γ.
- Paramètres NRTL/UNIQUAC de type ChemSep (à remplacer par ceux de ton cours si besoin) ; coûts = ordres de grandeur.