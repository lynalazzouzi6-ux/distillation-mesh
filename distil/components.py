"""
Base de données des corps purs (valeurs de manuel, ordre de grandeur) :
  antoine : (A, B, C)  log10(P[mmHg]) = A - B/(T[°C] + C)
  cp_l    : Cp liquide moyen (kJ/mol/K)
  lv_nb   : chaleur latente à la température normale d'ébullition (kJ/mol)
  tnb, tc : températures normale d'ébullition et critique (K)
  r, q    : paramètres de volume et de surface UNIQUAC
Tu peux ajouter un corps pur en ajoutant une entrée dans le dictionnaire.
"""

COMPOSANTS = {
    "ethanol": dict(nom="Éthanol", antoine=(8.20417, 1642.89, 230.300), cp_l=0.1250,
                    lv_nb=38.56, tnb=351.4, tc=513.9, r=2.1055, q=1.972),
    "eau":     dict(nom="Eau", antoine=(8.07131, 1730.63, 233.426), cp_l=0.0753,
                    lv_nb=40.66, tnb=373.15, tc=647.1, r=0.92, q=1.40),
    "methanol": dict(nom="Méthanol", antoine=(8.08097, 1582.271, 239.726), cp_l=0.0870,
                     lv_nb=35.21, tnb=337.7, tc=512.6, r=1.4311, q=1.432),
    "acetone": dict(nom="Acétone", antoine=(7.11714, 1210.595, 229.664), cp_l=0.1300,
                    lv_nb=29.10, tnb=329.2, tc=508.1, r=2.5735, q=2.336),
    "ipa":     dict(nom="Isopropanol", antoine=(8.11778, 1580.92, 219.61), cp_l=0.1700,
                    lv_nb=39.85, tnb=355.4, tc=508.3, r=3.2491, q=3.124),
}
