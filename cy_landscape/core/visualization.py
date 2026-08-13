"""
visualization.py — Visualisation du paysage CY.

Produit :
  1. Le Hodge shield (h11 vs h21) avec coloration par nombre de générations
  2. La distribution des générations
  3. Le Hodge shield filtré (candidats SM mis en évidence)
"""

import matplotlib
matplotlib.use("Agg")  # backend non-interactif
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
from typing import List, Optional
from cy_landscape.core.database import CYDatabase, CalabiYauManifold


# ─── Palette de couleurs ──────────────────────────────────────────

# Colormap personnalisée : bleu (n_gen éloigné) → rouge (n_gen = 3)
SM_CMAP = plt.cm.coolwarm_r


def plot_hodge_shield(
    db: CYDatabase,
    highlight_ngen: Optional[int] = 3,
    output_path: str = "hodge_shield.png",
    figsize: tuple = (14, 10),
    dpi: int = 150,
) -> str:
    """
    Trace le Hodge shield : nuage de points (h11, h21) coloré par n_gen.

    Le Hodge shield est la figure emblématique de la classification CY.
    Sa forme caractéristique en "bouclier" reflète les contraintes
    topologiques sur les nombres de Hodge.

    Args:
        db : base de données CY chargée
        highlight_ngen : nombre de générations à mettre en évidence (défaut: 3)
        output_path : chemin du fichier de sortie
        figsize : taille de la figure
        dpi : résolution

    Returns:
        Chemin du fichier sauvegardé.
    """
    fig, ax = plt.subplots(figsize=figsize, facecolor="#0a0a1a")
    ax.set_facecolor("#0a0a1a")

    # ── Données ──
    all_h11 = np.array([m.h11 for m in db.manifolds])
    all_h21 = np.array([m.h21 for m in db.manifolds])
    all_ngen = np.array([m.n_gen for m in db.manifolds])
    all_count = np.array([m.count for m in db.manifolds])

    # Taille des points proportionnelle au log du nombre de variétés
    sizes = 15 + 40 * np.log1p(all_count)

    # ── Couche 1 : tous les points, colorés par n_gen ──
    ngen_max = max(all_ngen) if len(all_ngen) > 0 else 1
    norm = mcolors.Normalize(vmin=0, vmax=ngen_max)

    scatter = ax.scatter(
        all_h11, all_h21,
        c=all_ngen, cmap="viridis", norm=norm,
        s=sizes, alpha=0.5, edgecolors="none",
        zorder=1,
    )

    # ── Couche 2 : mise en évidence des candidats SM ──
    if highlight_ngen is not None:
        mask = all_ngen == highlight_ngen
        if mask.any():
            ax.scatter(
                all_h11[mask], all_h21[mask],
                c="#ff3366", s=sizes[mask] * 2.5,
                alpha=0.95, edgecolors="white", linewidths=0.8,
                zorder=2, label=f"N_gen = {highlight_ngen} (Modèle Standard)",
            )

    # ── Diagonale miroir ──
    max_val = max(max(all_h11), max(all_h21)) + 5
    ax.plot([0, max_val], [0, max_val], "--", color="#ffffff",
            alpha=0.2, linewidth=0.8, label="Symétrie miroir (h¹¹ = h²¹)")

    # ── Décoration ──
    ax.set_xlabel("h¹·¹  (déformations de Kähler)", fontsize=13, color="white",
                   labelpad=10)
    ax.set_ylabel("h²·¹  (déformations complexes)", fontsize=13, color="white",
                   labelpad=10)
    ax.set_title("Hodge Shield — Paysage des variétés de Calabi-Yau (CICY)",
                  fontsize=16, color="white", pad=20, fontweight="bold")

    ax.tick_params(colors="white", labelsize=10)
    for spine in ax.spines.values():
        spine.set_color("#333355")

    # Colorbar
    cbar = plt.colorbar(scatter, ax=ax, shrink=0.7, pad=0.02)
    cbar.set_label("Nombre de générations |h¹¹ − h²¹|", color="white",
                    fontsize=11, labelpad=10)
    cbar.ax.yaxis.set_tick_params(color="white")
    cbar.ax.tick_params(labelcolor="white", labelsize=9)

    # Légende
    legend = ax.legend(loc="upper right", fontsize=10, framealpha=0.3,
                        facecolor="#1a1a2e", edgecolor="#333355",
                        labelcolor="white")

    # Annotations
    total = sum(m.count for m in db.manifolds)
    sm_candidates = [m for m in db.manifolds if m.n_gen == highlight_ngen]
    sm_total = sum(m.count for m in sm_candidates)

    info_text = (
        f"Base : {len(db.manifolds)} topologies, {total} variétés\n"
        f"Candidats SM (N_gen={highlight_ngen}) : "
        f"{len(sm_candidates)} topologies, {sm_total} variétés"
    )
    ax.text(0.02, 0.97, info_text, transform=ax.transAxes,
            fontsize=9, color="#aaaacc", verticalalignment="top",
            fontfamily="monospace",
            bbox=dict(boxstyle="round,pad=0.5", facecolor="#0a0a1a",
                      edgecolor="#333355", alpha=0.8))

    plt.tight_layout()
    fig.savefig(output_path, dpi=dpi, facecolor=fig.get_facecolor(),
                bbox_inches="tight")
    plt.close(fig)
    return output_path


def plot_generation_distribution(
    db: CYDatabase,
    output_path: str = "generation_distribution.png",
    figsize: tuple = (14, 6),
    dpi: int = 150,
) -> str:
    """
    Histogramme de la distribution du nombre de générations.

    Montre combien de topologies CY produisent 0, 1, 2, 3, ... N générations,
    mettant en évidence la rareté (ou non) des variétés à 3 générations.
    """
    dist = db.generation_distribution()
    ngens = sorted(dist.keys())
    topo_counts = [dist[n]["topologies"] for n in ngens]
    var_counts = [dist[n]["varietes"] for n in ngens]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize, facecolor="#0a0a1a")

    for ax in (ax1, ax2):
        ax.set_facecolor("#0a0a1a")
        ax.tick_params(colors="white", labelsize=9)
        for spine in ax.spines.values():
            spine.set_color("#333355")

    # ── Graphe 1 : par topologies distinctes ──
    colors1 = ["#ff3366" if n == 3 else "#4488cc" for n in ngens]
    bars1 = ax1.bar(ngens, topo_counts, color=colors1, alpha=0.85, edgecolor="none")
    ax1.set_xlabel("Nombre de générations |h¹¹ − h²¹|", color="white", fontsize=11)
    ax1.set_ylabel("Topologies distinctes", color="white", fontsize=11)
    ax1.set_title("Distribution par topologies", color="white", fontsize=13,
                   fontweight="bold")

    # Annoter la barre n=3
    if 3 in dist:
        idx3 = ngens.index(3)
        ax1.annotate(
            f"← {dist[3]['topologies']} topologies\n   à 3 générations",
            xy=(3, topo_counts[idx3]),
            xytext=(3 + 8, topo_counts[idx3] * 0.9),
            fontsize=9, color="#ff3366",
            arrowprops=dict(arrowstyle="->", color="#ff3366", lw=1.5),
        )

    # ── Graphe 2 : par nombre de variétés (avec multiplicité) ──
    colors2 = ["#ff3366" if n == 3 else "#44cc88" for n in ngens]
    bars2 = ax2.bar(ngens, var_counts, color=colors2, alpha=0.85, edgecolor="none")
    ax2.set_xlabel("Nombre de générations |h¹¹ − h²¹|", color="white", fontsize=11)
    ax2.set_ylabel("Nombre de variétés (avec multiplicité)", color="white", fontsize=11)
    ax2.set_title("Distribution par variétés", color="white", fontsize=13,
                   fontweight="bold")

    if 3 in dist:
        idx3 = ngens.index(3)
        ax2.annotate(
            f"← {dist[3]['varietes']} variétés\n   à 3 générations",
            xy=(3, var_counts[idx3]),
            xytext=(3 + 8, var_counts[idx3] * 0.9),
            fontsize=9, color="#ff3366",
            arrowprops=dict(arrowstyle="->", color="#ff3366", lw=1.5),
        )

    fig.suptitle("Combien de variétés CY reproduisent N générations de fermions ?",
                 fontsize=15, color="white", fontweight="bold", y=1.02)

    plt.tight_layout()
    fig.savefig(output_path, dpi=dpi, facecolor=fig.get_facecolor(),
                bbox_inches="tight")
    plt.close(fig)
    return output_path


def plot_sm_candidates_detail(
    candidates: List[CalabiYauManifold],
    output_path: str = "sm_candidates.png",
    figsize: tuple = (12, 8),
    dpi: int = 150,
) -> str:
    """
    Visualisation détaillée des candidats Modèle Standard (n_gen = 3).

    Montre chaque candidat dans le plan (h11, h21) avec des informations
    sur la caractéristique d'Euler et la multiplicité.
    """
    if not candidates:
        print("Aucun candidat SM trouvé.")
        return ""

    fig, ax = plt.subplots(figsize=figsize, facecolor="#0a0a1a")
    ax.set_facecolor("#0a0a1a")

    h11_vals = [m.h11 for m in candidates]
    h21_vals = [m.h21 for m in candidates]
    counts = [m.count for m in candidates]
    eulers = [m.euler for m in candidates]

    sizes = [30 + 60 * np.log1p(c) for c in counts]
    colors_by_euler = ["#ff6688" if e > 0 else "#66aaff" for e in eulers]

    scatter = ax.scatter(
        h11_vals, h21_vals,
        s=sizes, c=colors_by_euler,
        alpha=0.9, edgecolors="white", linewidths=0.5,
    )

    # Annoter chaque point
    for m in candidates:
        label = f"({m.h11},{m.h21})"
        offset = (5, 5) if m.h11 <= m.h21 else (5, -10)
        ax.annotate(
            label, (m.h11, m.h21), textcoords="offset points",
            xytext=offset, fontsize=7, color="#ddddee", alpha=0.8,
        )

    # Diagonale
    max_val = max(max(h11_vals), max(h21_vals)) + 3
    ax.plot([0, max_val], [0, max_val], "--", color="white", alpha=0.15, lw=0.8)

    # Légende personnalisée
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker="o", color="#0a0a1a", markerfacecolor="#ff6688",
               markersize=10, label="χ > 0  (h¹¹ > h²¹)"),
        Line2D([0], [0], marker="o", color="#0a0a1a", markerfacecolor="#66aaff",
               markersize=10, label="χ < 0  (h¹¹ < h²¹)"),
    ]
    ax.legend(handles=legend_elements, loc="upper left", fontsize=10,
              framealpha=0.3, facecolor="#1a1a2e", edgecolor="#333355",
              labelcolor="white")

    ax.set_xlabel("h¹·¹", fontsize=13, color="white")
    ax.set_ylabel("h²·¹", fontsize=13, color="white")
    ax.set_title(
        f"Candidats Modèle Standard : {len(candidates)} topologies à 3 générations",
        fontsize=14, color="white", fontweight="bold", pad=15,
    )

    ax.tick_params(colors="white", labelsize=10)
    for spine in ax.spines.values():
        spine.set_color("#333355")

    info = (
        f"Variétés totales : {sum(m.count for m in candidates)}\n"
        f"h¹¹ ∈ [{min(h11_vals)}, {max(h11_vals)}]\n"
        f"h²¹ ∈ [{min(h21_vals)}, {max(h21_vals)}]"
    )
    ax.text(0.98, 0.02, info, transform=ax.transAxes,
            fontsize=9, color="#aaaacc", verticalalignment="bottom",
            horizontalalignment="right", fontfamily="monospace",
            bbox=dict(boxstyle="round,pad=0.5", facecolor="#0a0a1a",
                      edgecolor="#333355", alpha=0.8))

    plt.tight_layout()
    fig.savefig(output_path, dpi=dpi, facecolor=fig.get_facecolor(),
                bbox_inches="tight")
    plt.close(fig)
    return output_path
