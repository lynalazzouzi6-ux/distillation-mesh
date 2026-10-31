"""
Systèmes binaires disponibles (paramètres type ChemSep, tau_ij = b_ij / T pour NRTL,
tau_ij = exp(b_ij / T) pour UNIQUAC ; T en K).
Pour ajouter un système : copie un bloc, change les corps purs (voir components.py),
les paramètres et les spécifications par défaut.

defauts : F (kmol/h), zF, xD, xW (spécifications de pureté), q (état de l'alimentation).
"""
from .thermo_models import Systeme

# Données expérimentales éthanol(1)/eau(2) à 760 mmHg, type tableau de manuel (Perry/Carey-Lewis).
# (x1, y1, T°C) -- à vérifier avec ta propre référence avant de les citer dans un rapport.
EXP_ETHANOL_EAU = [
    (0.0190, 0.1700, 95.5), (0.0721, 0.3891, 89.0), (0.0966, 0.4375, 86.7),
    (0.1238, 0.4704, 85.3), (0.1661, 0.5089, 84.1), (0.2337, 0.5445, 82.7),
    (0.2608, 0.5580, 82.3), (0.3273, 0.5826, 81.5), (0.3965, 0.6122, 80.7),
    (0.5079, 0.6564, 79.8), (0.5732, 0.6841, 79.3), (0.6763, 0.7385, 78.74),
    (0.7472, 0.7815, 78.41), (0.8943, 0.8943, 78.15),
]


def charger(cle):
    if cle not in SYSTEMES:
        raise KeyError(f"Système '{cle}' inconnu. Choix : {', '.join(SYSTEMES)}")
    return SYSTEMES[cle]()


def _ethanol_eau():
    return Systeme("ethanol_eau", "ethanol", "eau", "NRTL",
                   dict(b12=-29.1667, b21=624.868, alpha=0.2937),
                   defauts=dict(F=100.0, zF=0.20, xD=0.80, xW=0.01, q=1.0),
                   donnees_exp=EXP_ETHANOL_EAU,
                   description="Éthanol/eau, NRTL, azéotrope vers x = 0,88")

def _ethanol_eau_uniquac():
    return Systeme("ethanol_eau_uniquac", "ethanol", "eau", "UNIQUAC",
                   dict(b12=-87.46005814, b21=-55.28807596),
                   defauts=dict(F=100.0, zF=0.20, xD=0.80, xW=0.01, q=1.0),
                   donnees_exp=EXP_ETHANOL_EAU,
                   description="Éthanol/eau, UNIQUAC, azéotrope vers x = 0,88")

def _methanol_eau():
    return Systeme("methanol_eau", "methanol", "eau", "NRTL",
                   dict(b12=-95.13209283, b21=398.95345259, alpha=0.2999),
                   defauts=dict(F=100.0, zF=0.40, xD=0.95, xW=0.02, q=1.0),
                   description="Méthanol/eau, NRTL, sans azéotrope")

def _acetone_methanol():
    return Systeme("acetone_methanol", "acetone", "methanol", "NRTL",
                   dict(b12=59.42031348, b21=149.07536490, alpha=0.3003),
                   defauts=dict(F=100.0, zF=0.30, xD=0.70, xW=0.02, q=1.0),
                   description="Acétone/méthanol, NRTL, azéotrope à pression atmosphérique")

SYSTEMES = {
    "ethanol_eau": _ethanol_eau,
    "ethanol_eau_uniquac": _ethanol_eau_uniquac,
    "methanol_eau": _methanol_eau,
    "acetone_methanol": _acetone_methanol,
}
