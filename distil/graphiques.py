"""Toutes les figures du projet (matplotlib, enregistrées en PNG)."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def _sauver(fig, chemin):
    fig.tight_layout(); fig.savefig(chemin, dpi=150); plt.close(fig)


def fig_vle(sys, chemin, validation=None):
    xx = np.linspace(0.0005, 0.9995, 300)
    pts = [sys.point_bulle(x) for x in xx]
    T = [p[0] for p in pts]; y = [p[1] for p in pts]
    g = np.array([sys.gammas(x, t) for x, t in zip(xx, T)])
    fig, ax = plt.subplots(1, 3, figsize=(17, 5))
    ax[0].plot(xx, y, "b-", label=f"y(x) {sys.modele}"); ax[0].plot([0, 1], [0, 1], "k-", lw=.8)
    if validation:
        ax[0].plot(validation["x"], validation["y_exp"], "ro", ms=4, label="données exp.")
    if sys.x_az: ax[0].plot(sys.x_az, sys.x_az, "k*", ms=12, label=f"azéotrope x={sys.x_az:.3f}")
    ax[0].set(xlabel=f"x ({sys.nomA})", ylabel=f"y ({sys.nomA})", title="Diagramme x-y"); ax[0].legend(); ax[0].grid(alpha=.3)
    ax[1].plot(xx, T, "b-", label="T bulle"); ax[1].plot(y, T, "r-", label="T rosée")
    if validation:
        ax[1].plot(validation["x"], validation["T_exp"], "bo", ms=4)
        ax[1].plot(validation["y_exp"], validation["T_exp"], "ro", ms=4)
    ax[1].set(xlabel=f"x, y ({sys.nomA})", ylabel="T (°C)", title="Diagramme T-x-y"); ax[1].legend(); ax[1].grid(alpha=.3)
    ax[2].semilogy(xx, g[:, 0], label=f"γ {sys.nomA}"); ax[2].semilogy(xx, g[:, 1], label=f"γ {sys.nomB}")
    ax[2].set(xlabel=f"x ({sys.nomA})", ylabel="γ", title="Coefficients d'activité"); ax[2].legend(); ax[2].grid(alpha=.3, which="both")
    fig.suptitle(sys.description)
    _sauver(fig, chemin)


def fig_ponchon(sys, res, plateaux, ops, titre, chemin):
    xs = np.linspace(0, sys.x_max, 300)
    fig, ax = plt.subplots(figsize=(9, 7))
    ax.plot(xs, [sys.h_sat(x) for x in xs], "b-", label="Ébullition h(x)")
    ax.plot(xs, [sys.H_sat(y) for y in xs], "r-", label="Rosée H(y)")
    for (x, h, y, H) in plateaux: ax.plot([x, y], [h, H], "g-", lw=1)
    for (pl, pv, pole) in ops: ax.plot([pl[0], pv[0]], [pl[1], pv[1]], "k--", lw=.8)
    ax.plot([res["xD"], res["xW"]], [res["H_DELTA"], res["H_DELTA_P"]], "m:", lw=1.2)
    for lab, (px, py), c in [("Δ", (res["xD"], res["H_DELTA"]), "m"), ("Δ'", (res["xW"], res["H_DELTA_P"]), "m"),
                             ("F", (res["zF"], res["hF"]), "orange")]:
        ax.plot(px, py, "o", color=c); ax.annotate(lab, (px, py), textcoords="offset points", xytext=(6, 6))
    ax.set(xlabel=f"x, y ({sys.nomA})", ylabel="h, H (kJ/mol)", title=titre); ax.legend(); ax.grid(alpha=.3)
    _sauver(fig, chemin)


def fig_profils(sol, sys, titre, chemin):
    j = np.arange(1, sol.N + 1)
    fig, ax = plt.subplots(1, 3, figsize=(16, 5))
    ax[0].plot(j, sol.x, "bo-", ms=4, label="x (liquide)"); ax[0].plot(j, sol.y, "rs-", ms=4, label="y (vapeur)")
    ax[0].axvline(sol.NF, color="gray", ls="--", label=f"alimentation (étage {sol.NF})")
    ax[0].set(xlabel="étage (1 = tête, N = rebouilleur)", ylabel=f"fraction molaire {sys.nomA}", title="Compositions"); ax[0].legend(); ax[0].grid(alpha=.3)
    ax[1].plot(j, sol.T, "k^-", ms=4); ax[1].axvline(sol.NF, color="gray", ls="--")
    ax[1].set(xlabel="étage", ylabel="T (°C)", title="Profil de température"); ax[1].grid(alpha=.3)
    ax[2].step(j, sol.L, "b-", where="mid", label="L (liquide)"); ax[2].step(j, sol.V, "r-", where="mid", label="V (vapeur)")
    ax[2].axvline(sol.NF, color="gray", ls="--")
    ax[2].set(xlabel="étage", ylabel="kmol/h", title="Débits (non constants : le MESH n'impose pas le CMO)"); ax[2].legend(); ax[2].grid(alpha=.3)
    fig.suptitle(titre)
    _sauver(fig, chemin)


def fig_equivalence(eq, sys, chemin):
    sol = eq["sol"]; pl = eq["plateaux"]
    j = np.arange(1, sol.N + 1)
    xps = np.array([p[0] for p in pl])[:sol.N]
    yps = np.array([p[2] for p in pl])[:sol.N]
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    ax[0].plot(j, sol.x, "bo", ms=8, label="MESH x"); ax[0].plot(j[:len(xps)], xps, "k+", ms=10, label="Ponchon-Savarit x")
    ax[0].plot(j, sol.y, "rs", ms=8, mfc="none", label="MESH y"); ax[0].plot(j[:len(yps)], yps, "k+", ms=10, label="Ponchon-Savarit y")
    ax[0].set(xlabel="étage", ylabel="fraction molaire", title="Profils MESH vs Ponchon-Savarit"); ax[0].legend(); ax[0].grid(alpha=.3)
    ax[1].semilogy(j[:len(xps)], np.abs(sol.x[:len(xps)] - xps) + 1e-16, "ko-")
    ax[1].set(xlabel="étage", ylabel="|x_MESH - x_PS|", title=f"Écart max = {eq['ecart_x']:.1e}"); ax[1].grid(alpha=.3, which="both")
    _sauver(fig, chemin)


def fig_effet_R_q(courbes, chemin):
    """courbes : liste de (q, R_min, R[], N[], Qr[])."""
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    for q, rm, R, N, Qr in courbes:
        ax[0].plot(R, N, label=f"q={q} (Rmin={rm:.2f})"); ax[1].plot(R, Qr, label=f"q={q}")
    ax[0].set(xlabel="R", ylabel="Étages théoriques (PS)", title="N en fonction de R", ylim=(0, 60)); ax[0].legend(); ax[0].grid(alpha=.3)
    ax[1].set(xlabel="R", ylabel="Qr (kW)", title="Charge du rebouilleur"); ax[1].legend(); ax[1].grid(alpha=.3)
    _sauver(fig, chemin)


def fig_optimisation(designs, opt_idx, eco, chemin):
    N = [d["cout"]["N"] for d in designs]
    fig, ax = plt.subplots(2, 3, figsize=(17, 9))
    ax = ax.ravel()
    ax[0].plot(N, [d["cout"]["capex_annuel"] / 1e3 for d in designs], "o-", label="CAPEX annualisé")
    ax[0].plot(N, [d["cout"]["opex"] / 1e3 for d in designs], "s-", label="OPEX")
    ax[0].plot(N, [d["cout"]["tac"] / 1e3 for d in designs], "k-", lw=2, label="Total (TAC)")
    ax[0].axvline(N[opt_idx], color="r", ls="--", label=f"optimum N={N[opt_idx]}")
    ax[0].set(xlabel="N (étages théoriques)", ylabel="k$/an", title="Coûts annuels"); ax[0].legend(); ax[0].grid(alpha=.3)
    ax[1].plot(N, [d["cout"]["R"] for d in designs], "o-"); ax[1].axvline(N[opt_idx], color="r", ls="--")
    ax[1].set(xlabel="N", ylabel="R requis", title="Reflux requis pour la pureté"); ax[1].grid(alpha=.3)
    ax[2].plot(N, [d["cout"]["NF"] for d in designs], "o-"); ax[2].axvline(N[opt_idx], color="r", ls="--")
    ax[2].set(xlabel="N", ylabel="étage d'alimentation optimal", title="Position de l'alimentation"); ax[2].grid(alpha=.3)
    cap = np.array([d["cout"]["capex_annuel"] / 1e3 for d in designs]); op = np.array([d["cout"]["opex"] / 1e3 for d in designs])
    sc = ax[3].scatter(cap, op, c=N, cmap="viridis", s=40)
    ax[3].plot(cap[opt_idx], op[opt_idx], "r*", ms=16, label="optimum TAC")
    ax[3].set(xlabel="CAPEX annualisé (k$/an)", ylabel="OPEX (k$/an)", title="Compromis CAPEX / OPEX (couleur = N)"); ax[3].legend(); ax[3].grid(alpha=.3)
    fig.colorbar(sc, ax=ax[3], label="N")
    ax[4].plot(N, [d["cout"]["Qr"] for d in designs], "o-", label="Qr"); ax[4].plot(N, [d["cout"]["Qc"] for d in designs], "s-", label="Qc")
    ax[4].set(xlabel="N", ylabel="kW", title="Charges thermiques"); ax[4].legend(); ax[4].grid(alpha=.3)
    ax[5].plot(N, [d["cout"]["diam"] for d in designs], "o-", label="diamètre (m)")
    ax2 = ax[5].twinx(); ax2.plot(N, [d["cout"]["haut"] for d in designs], "rs-", label="hauteur (m)")
    ax[5].set(xlabel="N", ylabel="diamètre (m)", title="Dimensions de la colonne"); ax2.set_ylabel("hauteur (m)"); ax[5].grid(alpha=.3)
    _sauver(fig, chemin)


def fig_sensibilite(prix, ropt_sur_rmin, nopt, tac, chemin):
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))
    ax[0].semilogx(prix, nopt, "go-"); ax[0].set(xlabel="Prix de la vapeur ($/t)", ylabel="N optimal", title="N optimal")
    ax[1].semilogx(prix, ropt_sur_rmin, "ro-"); ax[1].set(xlabel="Prix de la vapeur ($/t)", ylabel="R_opt / R_min", title="Reflux optimal")
    ax[2].semilogx(prix, np.array(tac) / 1e3, "ko-"); ax[2].set(xlabel="Prix de la vapeur ($/t)", ylabel="TAC (k$/an)", title="Coût total minimal")
    for a in ax: a.grid(alpha=.3, which="both")
    _sauver(fig, chemin)
