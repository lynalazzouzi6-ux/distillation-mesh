# Distillation Acétone/Méthanol — Acétone/méthanol, NRTL, azéotrope à pression atmosphérique
Modèle : NRTL | P = 1.00 atm | F = 100.0 kmol/h, zF = 0.3, q = 1.0
Spécifications : xD ≥ 0.7, xW ≤ 0.02
Azéotrope du modèle : x = 0.7903, T = 55.39 °C  (xD doit rester en dessous)

## 1. Validation de l'équilibre liquide-vapeur
Pas de données expérimentales enregistrées pour ce système (ajoute-les dans systemes.py).
Comparaison avec la bibliothèque thermo (x=0.3, 70 °C) : γ nous = 1.32205/1.05710, thermo = 1.32205/1.05710 (écart max 0.0e+00) ; hE nous = 0.3355, thermo = 0.3355 kJ/mol

## 2. Ponchon-Savarit
R_min = 2.296 | à R = 1,5 R_min = 3.444 : 16.63 étages théoriques, alimentation étage 10, Qc = 1568 kW, Qr = 1568 kW
  q = 1.0: R_min = 2.296
  q = 0.5: R_min = 2.800
  q = 0.0: R_min = 3.563

## 3. MESH rigoureux et équivalence avec Ponchon-Savarit
MESH à N = 17, N_F = 10, R = 3.444 : xD = 0.70207, xW = 0.01855, convergé en 4 itérations (résidu 2.6e-13)
Ponchon-Savarit refait avec ces puretés : 17.000 étages (attendu 17) ; écart max sur x = 8.0e-09, écart relatif Qc = 2.1e-10, Qr = 1.4e-10
Fermeture des bilans : matière 1.9e-13, composant 3.1e-15, énergie (relative) 4.8e-14

## 4. Optimisation CAPEX + OPEX (MESH, variables N, N_F, R)
31 conceptions évaluées (560 résolutions MESH).
**Optimum : N = 26 étages théoriques (36 plateaux réels, E_o = 0.7), alimentation étage 16, R = 2.489 (1.08 × R_min)**
Colonne : Ø 1.31 m, H 24.6 m | Qc = 1231 kW, Qr = 1231 kW | xD = 0.7000, xW = 0.0200
CAPEX = 1073 k$ (175 k$/an) | OPEX = 458 k$/an | TAC = 633 k$/an
Comparaison au même R : Ponchon-Savarit donne 26.00 étages (MESH : 26 entiers).
Vérification NLP (SLSQP sur R et D, contraintes xD ≥ 0.7, xW ≤ 0.02) : R = 2.4888, D = 41.176 kmol/h, TAC = 633.0 k$/an, succès = True -> optimum confirmé

## 5. Sensibilité au prix de la vapeur (mêmes conceptions, coûts recalculés)
  vapeur    5 $/t -> N_opt = 21, N_F = 13, R_opt = 2.77 (1.20 × R_min), TAC = 290 k$/an
  vapeur   10 $/t -> N_opt = 24, N_F = 15, R_opt = 2.57 (1.12 × R_min), TAC = 379 k$/an
  vapeur   25 $/t -> N_opt = 26, N_F = 16, R_opt = 2.49 (1.08 × R_min), TAC = 633 k$/an
  vapeur   50 $/t -> N_opt = 29, N_F = 18, R_opt = 2.41 (1.05 × R_min), TAC = 1042 k$/an
  vapeur  100 $/t -> N_opt = 33, N_F = 20, R_opt = 2.35 (1.03 × R_min), TAC = 1849 k$/an
  vapeur  200 $/t -> N_opt = 36, N_F = 22, R_opt = 2.33 (1.02 × R_min), TAC = 3440 k$/an

Temps de calcul : 11 s

## Limites du modèle
- Étages théoriques, P constante (pas de perte de charge), condenseur total, une alimentation, mélange binaire.
- Cp liquides constants, chaleurs latentes par la corrélation de Watson, vapeur idéale, enthalpie d'excès via le modèle de γ.
- Paramètres NRTL/UNIQUAC de type ChemSep (à remplacer par ceux de ton cours si besoin) ; coûts = ordres de grandeur.