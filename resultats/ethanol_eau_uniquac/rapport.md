# Distillation Éthanol/Eau — Éthanol/eau, UNIQUAC, azéotrope vers x = 0,88
Modèle : UNIQUAC | P = 1.00 atm | F = 100.0 kmol/h, zF = 0.2, q = 1.0
Spécifications : xD ≥ 0.8, xW ≤ 0.01
Azéotrope du modèle : x = 0.8782, T = 78.07 °C  (xD doit rester en dessous)

## 1. Validation de l'équilibre liquide-vapeur
Écart avec 14 points expérimentaux (données de manuel) : RMSE(y) = 0.0094, RMSE(T) = 0.23 °C
Comparaison avec la bibliothèque thermo (x=0.3, 70 °C) : γ nous = 1.75415/1.18925, thermo = 1.75415/1.18925 (écart max 6.7e-16) ; hE nous = 0.3863, thermo = 0.3863 kJ/mol

## 2. Ponchon-Savarit
R_min = 0.917 | à R = 1,5 R_min = 1.375 : 14.72 étages théoriques, alimentation étage 12, Qc = 618 kW, Qr = 631 kW
  q = 1.0: R_min = 0.917
  q = 0.5: R_min = 1.703
  q = 0.0: R_min = 3.592

## 3. MESH rigoureux et équivalence avec Ponchon-Savarit
MESH à N = 15, N_F = 12, R = 1.375 : xD = 0.80068, xW = 0.00978, convergé en 4 itérations (résidu 7.3e-11)
Ponchon-Savarit refait avec ces puretés : 15.000 étages (attendu 15) ; écart max sur x = 2.7e-08, écart relatif Qc = 3.2e-10, Qr = 7.5e-10
Fermeture des bilans : matière 4.3e-14, composant 6.7e-14, énergie (relative) 3.0e-14

## 4. Optimisation CAPEX + OPEX (MESH, variables N, N_F, R)
31 conceptions évaluées (531 résolutions MESH).
**Optimum : N = 19 étages théoriques (26 plateaux réels, E_o = 0.7), alimentation étage 16, R = 1.175 (1.28 × R_min)**
Colonne : Ø 0.82 m, H 18.6 m | Qc = 566 kW, Qr = 579 kW | xD = 0.8000, xW = 0.0100
CAPEX = 482 k$ (78 k$/an) | OPEX = 215 k$/an | TAC = 294 k$/an
Comparaison au même R : Ponchon-Savarit donne 19.00 étages (MESH : 19 entiers).
Vérification NLP (SLSQP sur R et D, contraintes xD ≥ 0.8, xW ≤ 0.01) : R = 1.1752, D = 24.051 kmol/h, TAC = 293.7 k$/an, succès = True -> optimum confirmé

## 5. Sensibilité au prix de la vapeur (mêmes conceptions, coûts recalculés)
  vapeur    5 $/t -> N_opt = 15, N_F = 13, R_opt = 1.34 (1.47 × R_min), TAC = 131 k$/an
  vapeur   10 $/t -> N_opt = 15, N_F = 13, R_opt = 1.34 (1.47 × R_min), TAC = 173 k$/an
  vapeur   25 $/t -> N_opt = 19, N_F = 16, R_opt = 1.18 (1.28 × R_min), TAC = 294 k$/an
  vapeur   50 $/t -> N_opt = 24, N_F = 21, R_opt = 1.08 (1.17 × R_min), TAC = 484 k$/an
  vapeur  100 $/t -> N_opt = 29, N_F = 26, R_opt = 1.03 (1.12 × R_min), TAC = 854 k$/an
  vapeur  200 $/t -> N_opt = 36, N_F = 33, R_opt = 0.99 (1.08 × R_min), TAC = 1577 k$/an

Temps de calcul : 14 s

## Limites du modèle
- Étages théoriques, P constante (pas de perte de charge), condenseur total, une alimentation, mélange binaire.
- Cp liquides constants, chaleurs latentes par la corrélation de Watson, vapeur idéale, enthalpie d'excès via le modèle de γ.
- Paramètres NRTL/UNIQUAC de type ChemSep (à remplacer par ceux de ton cours si besoin) ; coûts = ordres de grandeur.