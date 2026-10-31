"""
Modèle de coût d'ordre de grandeur (CAPEX + OPEX) d'une colonne de distillation.
Les coefficients sont des valeurs de départ raisonnables : REMPLACE-LES par ceux de ton cours
(Turton, Peters & Timmerhaus, indices CEPCI...) pour un rapport final.
"""
import math
from dataclasses import dataclass


@dataclass
class ParamsEco:
    E_o: float = 0.70               # efficacité globale de plateau
    pas_plateau: float = 0.6        # m
    h_extra: float = 3.0            # m (fond + tête)
    u_vap: float = 0.8              # m/s, vitesse superficielle (dimensionnement du diamètre)
    c_coque: float = 3000.0         # $/m² de virole (installé)
    c_plateau: float = 1500.0       # $/m² de section, par plateau
    c0_hx: float = 12000.0          # $ pour un échangeur de 10 m² (loi en (A/10)^0.65)
    u_cond: float = 0.8             # kW/m²/K
    u_reb: float = 1.2              # kW/m²/K
    t_eau_in: float = 25.0          # °C
    t_eau_out: float = 35.0         # °C
    t_vapeur: float = 144.0         # °C (vapeur BP ~ 4 bar)
    lv_vapeur: float = 2130.0       # kJ/kg
    prix_vapeur: float = 25.0       # $/t
    prix_eau: float = 0.05          # $/t
    heures: float = 8000.0          # h/an
    facteur_install: float = 2.5    # facteur d'installation (Lang simplifié)
    taux: float = 0.10              # taux d'actualisation
    duree: int = 10                 # ans

    @property
    def crf(self):
        i, n = self.taux, self.duree
        return i * (1 + i) ** n / ((1 + i) ** n - 1)


def evaluer(sol, eco: ParamsEco):
    """Coûts d'une solution MESH convergée. Retourne un dictionnaire."""
    N = sol.N
    n_reels = math.ceil(max(N - 1, 1) / eco.E_o)             # le rebouilleur compte pour 1 étage
    T_top, T_bot = float(sol.T[0]), float(sol.T[-1])
    # débit volumique de vapeur maximal (gaz parfait, P ~ 1 atm)
    q_vol = float(max(sol.V * (sol.T + 273.15)) * 22.414 / 273.15 / 3600.0)   # m³/s
    diam = math.sqrt(4 * q_vol / (math.pi * eco.u_vap))
    haut = n_reels * eco.pas_plateau + eco.h_extra
    c_col = (eco.c_coque * math.pi * diam * haut
             + eco.c_plateau * (math.pi * diam ** 2 / 4) * n_reels)
    Qc, Qr = sol.Qc / 3.6, sol.Qr / 3.6                       # MJ/h -> kW
    d1, d2 = T_top - eco.t_eau_out, T_top - eco.t_eau_in
    dtlm = (d1 - d2) / math.log(d1 / d2) if abs(d1 - d2) > 1e-9 else d1
    A_c = Qc / (eco.u_cond * dtlm)
    A_r = Qr / (eco.u_reb * (eco.t_vapeur - T_bot))
    c_hx = eco.c0_hx * ((A_c / 10) ** 0.65 + (A_r / 10) ** 0.65)
    capex = eco.facteur_install * (c_col + c_hx)
    vapeur_t_h = Qr / eco.lv_vapeur * 3.6
    eau_t_h = Qc / (4.18 * (eco.t_eau_out - eco.t_eau_in)) * 3.6
    opex = eco.heures * (vapeur_t_h * eco.prix_vapeur + eau_t_h * eco.prix_eau)
    return dict(N=N, NF=sol.NF, R=sol.R, n_reels=n_reels, diam=diam, haut=haut, Qc=Qc, Qr=Qr,
                A_c=A_c, A_r=A_r, capex=capex, capex_annuel=capex * eco.crf, opex=opex,
                tac=opex + capex * eco.crf)
