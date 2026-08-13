"""
viz_stage2.py — Visualisations pour l'étage 2 (fibrés vectoriels).
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from typing import List, Dict
from cy_landscape.core.bundles import BundleAnalysis


def plot_gauge_group_distribution(
    results: List[BundleAnalysis],
    output_path: str = "gauge_groups.png",
    figsize: tuple = (14, 6),
    dpi: int = 150,
) -> str:
    """
    Distribution des groupes de jauge parmi les fibrés viables.
    """
    # Compter par groupe
    group_counts = {}
    group_scores = {}
    for r in results:
        g = r.gauge_group
        group_counts[g] = group_counts.get(g, 0) + 1
        if g not in group_scores:
            group_scores[g] = []
        group_scores[g].append(r.sm_score)

    if not group_counts:
        return ""

    groups = sorted(group_counts.keys(), key=lambda g: group_counts[g], reverse=True)
    counts = [group_counts[g] for g in groups]
    avg_scores = [np.mean(group_scores[g]) for g in groups]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize, facecolor="#0a0a1a")

    for ax in (ax1, ax2):
        ax.set_facecolor("#0a0a1a")
        ax.tick_params(colors="white", labelsize=10)
        for spine in ax.spines.values():
            spine.set_color("#333355")

    # Couleurs par viabilité GUT
    gut_colors = {
        "E₆": "#ff6b6b", "SO(10)": "#ffd93d", "SU(5)": "#6bcb77",
        "E₇": "#4d96ff", "E₈": "#8b8b8b",
    }
    colors = [gut_colors.get(g, "#aa88cc") for g in groups]

    # ── Graphe 1 : nombre de fibrés par groupe ──
    bars = ax1.barh(range(len(groups)), counts, color=colors, alpha=0.85)
    ax1.set_yticks(range(len(groups)))
    ax1.set_yticklabels(groups, fontsize=12, color="white")
    ax1.set_xlabel("Nombre de fibrés viables", color="white", fontsize=11)
    ax1.set_title("Fibrés par groupe de jauge", color="white",
                   fontsize=13, fontweight="bold")
    ax1.invert_yaxis()

    for i, (cnt, g) in enumerate(zip(counts, groups)):
        ax1.text(cnt + 0.5, i, str(cnt), va="center", color="white", fontsize=10)

    # ── Graphe 2 : score SM moyen par groupe ──
    bars2 = ax2.barh(range(len(groups)), avg_scores, color=colors, alpha=0.85)
    ax2.set_yticks(range(len(groups)))
    ax2.set_yticklabels(groups, fontsize=12, color="white")
    ax2.set_xlabel("Score SM moyen", color="white", fontsize=11)
    ax2.set_title("Compatibilité Modèle Standard", color="white",
                   fontsize=13, fontweight="bold")
    ax2.set_xlim(0, 105)
    ax2.invert_yaxis()

    for i, sc in enumerate(avg_scores):
        ax2.text(sc + 1, i, f"{sc:.0f}", va="center", color="white", fontsize=10)

    # Légende
    legend_patches = [
        mpatches.Patch(color="#6bcb77", label="SU(5) — le plus proche du SM"),
        mpatches.Patch(color="#ffd93d", label="SO(10) — GUT classique"),
        mpatches.Patch(color="#ff6b6b", label="E₆ — GUT classique"),
    ]
    fig.legend(handles=legend_patches, loc="lower center", ncol=3,
               fontsize=10, framealpha=0.3, facecolor="#1a1a2e",
               edgecolor="#333355", labelcolor="white",
               bbox_to_anchor=(0.5, -0.02))

    plt.tight_layout()
    fig.savefig(output_path, dpi=dpi, facecolor=fig.get_facecolor(),
                bbox_inches="tight")
    plt.close(fig)
    return output_path


def plot_score_landscape(
    results: List[BundleAnalysis],
    output_path: str = "score_landscape.png",
    figsize: tuple = (14, 8),
    dpi: int = 150,
) -> str:
    """
    Paysage des scores SM dans le plan (h11, h21), coloré par meilleur score.
    """
    if not results:
        return ""

    # Meilleur score par (h11, h21)
    best_by_cy = {}
    for r in results:
        key = (r.geometry.h11, r.geometry.h21)
        if key not in best_by_cy or r.sm_score > best_by_cy[key].sm_score:
            best_by_cy[key] = r

    fig, ax = plt.subplots(figsize=figsize, facecolor="#0a0a1a")
    ax.set_facecolor("#0a0a1a")

    h11_vals = [k[0] for k in best_by_cy]
    h21_vals = [k[1] for k in best_by_cy]
    scores = [best_by_cy[k].sm_score for k in best_by_cy]
    gauge_groups = [best_by_cy[k].gauge_group for k in best_by_cy]

    # Taille proportionnelle au score
    sizes = [30 + s * 3 for s in scores]

    scatter = ax.scatter(
        h11_vals, h21_vals,
        c=scores, cmap="RdYlGn", vmin=0, vmax=100,
        s=sizes, alpha=0.9, edgecolors="white", linewidths=0.5,
    )

    # Annoter avec le groupe de jauge
    for key in best_by_cy:
        r = best_by_cy[key]
        label = f"{r.gauge_group}\n{r.sm_score:.0f}"
        ax.annotate(
            label, (key[0], key[1]),
            textcoords="offset points", xytext=(8, -5),
            fontsize=7, color="#ddddee", alpha=0.8,
        )

    # Diagonale miroir
    max_val = max(max(h11_vals), max(h21_vals)) + 3
    ax.plot([0, max_val], [0, max_val], "--", color="white", alpha=0.15, lw=0.8)

    cbar = plt.colorbar(scatter, ax=ax, shrink=0.7)
    cbar.set_label("Score SM", color="white", fontsize=11)
    cbar.ax.tick_params(labelcolor="white")

    ax.set_xlabel("h¹·¹", fontsize=13, color="white")
    ax.set_ylabel("h²·¹", fontsize=13, color="white")
    ax.set_title(
        "Paysage des scores — Meilleur fibré par topologie CY",
        fontsize=15, color="white", fontweight="bold", pad=15,
    )
    ax.tick_params(colors="white", labelsize=10)
    for spine in ax.spines.values():
        spine.set_color("#333355")

    plt.tight_layout()
    fig.savefig(output_path, dpi=dpi, facecolor=fig.get_facecolor(),
                bbox_inches="tight")
    plt.close(fig)
    return output_path


def plot_bundle_summary_table(
    top_results: List[BundleAnalysis],
    output_path: str = "top_bundles.png",
    figsize: tuple = (16, 10),
    dpi: int = 150,
) -> str:
    """
    Tableau visuel des meilleurs fibrés candidats.
    """
    if not top_results:
        return ""

    n = min(len(top_results), 15)
    results = top_results[:n]

    fig, ax = plt.subplots(figsize=figsize, facecolor="#0a0a1a")
    ax.set_facecolor("#0a0a1a")
    ax.axis("off")

    # En-têtes
    headers = ["#", "CY (h¹¹,h²¹)", "χ", "Rang", "Groupe", "c₁=0",
               "N_gen", "Anomalie", "Score"]
    col_widths = [0.04, 0.12, 0.06, 0.06, 0.10, 0.06, 0.08, 0.08, 0.08]
    col_x = [sum(col_widths[:i]) + 0.15 for i in range(len(headers))]

    # Dessiner les en-têtes
    y_start = 0.95
    dy = 0.055
    for i, (header, x) in enumerate(zip(headers, col_x)):
        ax.text(x, y_start, header, fontsize=10, fontweight="bold",
                color="#88aaff", transform=ax.transAxes, ha="center")

    # Ligne de séparation
    ax.plot([0.1, 0.9], [y_start - 0.015, y_start - 0.015],
            color="#333355", linewidth=0.8, transform=ax.transAxes)

    # Dessiner les lignes
    for row, r in enumerate(results):
        y = y_start - (row + 1) * dy
        check = lambda b: "✓" if b else "✗"
        color_check = lambda b: "#6bcb77" if b else "#ff6b6b"

        values = [
            (str(row + 1), "white"),
            (f"({r.geometry.h11},{r.geometry.h21})", "white"),
            (str(r.geometry.euler), "#aaaacc"),
            (f"SU({r.structure_group_rank})", "#ffcc66"),
            (r.gauge_group, "#88ddff"),
            (check(r.c1_vanishes), color_check(r.c1_vanishes)),
            (f"{r.chiral_index_value:.1f}", color_check(r.correct_generations)),
            (check(r.anomaly_cancelled), color_check(r.anomaly_cancelled)),
            (f"{r.sm_score:.0f}", "#6bcb77" if r.sm_score >= 70 else
             "#ffd93d" if r.sm_score >= 50 else "#ff6b6b"),
        ]

        # Fond alterné
        if row % 2 == 0:
            rect_bg = plt.Rectangle(
                (0.1, y - dy/2 + 0.01), 0.8, dy - 0.015,
                facecolor="#111133", alpha=0.5, transform=ax.transAxes,
            )
            ax.add_patch(rect_bg)

        for (val, color), x in zip(values, col_x):
            ax.text(x, y, val, fontsize=10, color=color,
                    transform=ax.transAxes, ha="center", va="center")

    ax.set_title(
        f"Top {n} fibrés vectoriels — Candidats Modèle Standard",
        fontsize=14, color="white", fontweight="bold", pad=20,
    )

    fig.savefig(output_path, dpi=dpi, facecolor=fig.get_facecolor(),
                bbox_inches="tight")
    plt.close(fig)
    return output_path


def plot_physics_chain(
    best_result: BundleAnalysis,
    output_path: str = "physics_chain.png",
    figsize: tuple = (16, 6),
    dpi: int = 150,
) -> str:
    """
    Schéma de la chaîne physique : CY → Fibré → Groupe de jauge → Particules.
    Illustration du meilleur candidat trouvé.
    """
    if best_result is None:
        return ""

    fig, ax = plt.subplots(figsize=figsize, facecolor="#0a0a1a")
    ax.set_facecolor("#0a0a1a")
    ax.axis("off")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    # ── Boîtes de la chaîne ──
    boxes = [
        {
            "title": "Variété CY",
            "content": (
                f"h¹¹ = {best_result.geometry.h11}\n"
                f"h²¹ = {best_result.geometry.h21}\n"
                f"χ = {best_result.geometry.euler}"
            ),
            "color": "#2a4a7f",
            "x": 0.08,
        },
        {
            "title": "Fibré vectoriel",
            "content": (
                f"Rang {best_result.structure_group_rank}\n"
                f"Structure : SU({best_result.structure_group_rank})\n"
                f"c₁ = 0 : {'✓' if best_result.c1_vanishes else '✗'}"
            ),
            "color": "#4a2a7f",
            "x": 0.30,
        },
        {
            "title": "Groupe de jauge",
            "content": (
                f"{best_result.gauge_group}\n"
                f"GUT viable : {'✓' if best_result.gut_viable else '✗'}\n"
                f"(commutant dans E₈)"
            ),
            "color": "#7f4a2a",
            "x": 0.52,
        },
        {
            "title": "Physique",
            "content": (
                f"N_gen = {best_result.chiral_index_value:.1f}\n"
                f"3 familles : {'✓' if best_result.correct_generations else '✗'}\n"
                f"Score SM : {best_result.sm_score:.0f}/100"
            ),
            "color": "#2a7f4a",
            "x": 0.74,
        },
    ]

    box_w, box_h = 0.18, 0.55
    box_y = 0.25

    for box in boxes:
        # Rectangle
        rect = plt.Rectangle(
            (box["x"], box_y), box_w, box_h,
            facecolor=box["color"], edgecolor="white",
            linewidth=1.5, alpha=0.8, transform=ax.transAxes,
        )
        ax.add_patch(rect)

        # Titre
        ax.text(
            box["x"] + box_w / 2, box_y + box_h - 0.06,
            box["title"], fontsize=11, fontweight="bold",
            color="white", ha="center", va="top", transform=ax.transAxes,
        )

        # Contenu
        ax.text(
            box["x"] + box_w / 2, box_y + box_h / 2 - 0.05,
            box["content"], fontsize=10,
            color="#ddddee", ha="center", va="center",
            transform=ax.transAxes, linespacing=1.6,
        )

    # Flèches entre les boîtes
    arrow_props = dict(
        arrowstyle="->", color="#ffcc66", lw=2,
        connectionstyle="arc3,rad=0",
    )
    for i in range(len(boxes) - 1):
        x_start = boxes[i]["x"] + box_w
        x_end = boxes[i + 1]["x"]
        y_mid = box_y + box_h / 2
        ax.annotate(
            "", xy=(x_end, y_mid), xytext=(x_start, y_mid),
            arrowprops=arrow_props, transform=ax.transAxes,
        )

    # Labels des flèches
    arrow_labels = ["compactification", "brisure de E₈", "spectre de matière"]
    for i, label in enumerate(arrow_labels):
        x_mid = (boxes[i]["x"] + box_w + boxes[i + 1]["x"]) / 2
        ax.text(
            x_mid, box_y + box_h / 2 + 0.08, label,
            fontsize=8, color="#ffcc66", ha="center", va="bottom",
            transform=ax.transAxes, fontstyle="italic",
        )

    ax.set_title(
        "Chaîne physique : de la géométrie aux particules",
        fontsize=15, color="white", fontweight="bold",
        y=0.95,
    )

    fig.savefig(output_path, dpi=dpi, facecolor=fig.get_facecolor(),
                bbox_inches="tight")
    plt.close(fig)
    return output_path
