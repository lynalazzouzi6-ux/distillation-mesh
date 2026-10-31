# Distillation Méthanol/Eau — Méthanol/eau, NRTL, sans azéotrope
Modèle : NRTL | P = 1.00 atm | F = 100.0 kmol/h, zF = 0.4, q = 1.0
Spécifications : xD ≥ 0.95, xW ≤ 0.02
Pas d'azéotrope : séparation possible jusqu'à la pureté complète.

## 1. Validation de l'équilibre liquide-vapeur
Pas de données expérimentales enregistrées pour ce système (ajoute-les dans systemes.py).
Comparaison avec la bibliothèque thermo (x=0.3, 70 °C) : γ nous = 1.36311/1.09324, thermo = 1.36311/1.09324 (écart max 0.0e+00) ; hE nous = 0.3516, thermo = 0.3516 kJ/mol

## 2. Ponchon-Savarit
R_min = 0.699 | à R = 1,5 R_min = 1.048 : 8.97 étages théoriques, alimentation étage 6, Qc = 825 kW, Qr = 833 kW
  q = 1.0: R_min = 0.699
  q = 0.5: R_min = 1.062
  q = 0.0: R_min = 2.023

## 3. MESH rigoureux et équivalence avec Ponchon-Savarit
MESH à N = 9, N_F = 6, R = 1.048 : xD = 0.95061, xW = 0.01958, convergé en 4 itérations (résidu 6.3e-13)
Ponchon-Savarit refait avec ces puretés : 9.000 étages (attendu 9) ; écart max sur x = 2.3e-08, écart relatif Qc = 7.9e-13, Qr = 7.6e-12
Fermeture des bilans : matière 1.4e-14, composant 1.3e-12, énergie (relative) 1.3e-14

## 4. Optimisation CAPEX + OPEX (MESH, variables N, N_F, R)
31 conceptions évaluées (522 résolutions MESH).
**Optimum : N = 15 étages théoriques (20 plateaux réels, E_o = 0.7), alimentation étage 10, R = 0.748 (1.07 × R_min)**
Colonne : Ø 0.94 m, H 15.0 m | Qc = 704 kW, Qr = 712 kW | xD = 0.9500, xW = 0.0200
CAPEX = 471 k$ (77 k$/an) | OPEX = 265 k$/an | TAC = 342 k$/an
Comparaison au même R : Ponchon-Savarit donne 15.00 étages (MESH : 15 entiers).
Vérification NLP (SLSQP sur R et D, contraintes xD ≥ 0.95, xW ≤ 0.02) : R = 0.7477, D = 40.860 kmol/h, TAC = 341.7 k$/an, succès = True -> optimum confirmé

## 5. Sensibilité au prix de la vapeur (mêmes conceptions, coûts recalculés)
  vapeur    5 $/t -> N_opt = 10, N_F = 7, R_opt = 0.94 (1.34 × R_min), TAC = 142 k$/an
  vapeur   10 $/t -> N_opt = 12, N_F = 8, R_opt = 0.82 (1.17 × R_min), TAC = 193 k$/an
  vapeur   25 $/t -> N_opt = 15, N_F = 10, R_opt = 0.75 (1.07 × R_min), TAC = 342 k$/an
  vapeur   50 $/t -> N_opt = 15, N_F = 10, R_opt = 0.75 (1.07 × R_min), TAC = 582 k$/an
  vapeur  100 $/t -> N_opt = 17, N_F = 12, R_opt = 0.73 (1.04 × R_min), TAC = 1059 k$/an
  vapeur  200 $/t -> N_opt = 21, N_F = 15, R_opt = 0.71 (1.01 × R_min), TAC = 2004 k$/an

Temps de calcul : 10 s

## Limites du modèle
- Étages théoriques, P constante (pas de perte de charge), condenseur total, une alimentation, mélange binaire.
- Cp liquides constants, chaleurs latentes par la corrélation de Watson, vapeur idéale, enthalpie d'excès via le modèle de γ.
- Paramètres NRTL/UNIQUAC de type ChemSep (à remplacer par ceux de ton cours si besoin) ; coûts = ordres de grandeur.