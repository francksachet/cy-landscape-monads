"""viz_stage3.py -- Visualisations pour l'etage 3."""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from typing import List


def plot_particle_spectrum(analyses, output_path="particle_spectra.png",
                           figsize=(16, 10), dpi=150):
    if not analyses:
        return ""
    n = min(len(analyses), 6)
    top = analyses[:n]
    rows = 2 if n > 3 else 1
    cols = min(n, 3)
    fig, axes = plt.subplots(rows, cols, figsize=figsize, facecolor="#0a0a1a")
    if rows == 1 and cols == 1:
        axes = np.array([axes])
    axes = axes.flat
    fig.suptitle("Spectres de particules - Top candidats vs Modele Standard",
                 fontsize=15, color="white", fontweight="bold", y=1.02)

    color_map = {"10": "#ff6b6b", "10bar": "#ff9999",
                 "5bar": "#ffd93d", "5": "#ffed8a",
                 "16": "#ff6b6b", "16bar": "#ff9999",
                 "27": "#ff6b6b", "27bar": "#ff9999",
                 "1": "#8b8b8b"}

    for idx in range(rows * cols):
        ax = axes[idx]
        ax.set_facecolor("#0a0a1a")
        ax.tick_params(colors="white", labelsize=8)
        for spine in ax.spines.values():
            spine.set_color("#333355")

        if idx >= n:
            ax.set_visible(False)
            continue

        analysis = top[idx]
        sp = analysis.spectrum
        ba = analysis.bundle_analysis
        reps = {k: v for k, v in sp.representations.items() if v > 0}
        if not reps:
            ax.text(0.5, 0.5, "Spectre vide", ha="center", va="center",
                    color="#666688", fontsize=10, transform=ax.transAxes)
            continue

        labels = list(reps.keys())
        values = list(reps.values())
        colors = [color_map.get(l, "#6bcb77") for l in labels]
        ax.bar(range(len(labels)), values, color=colors, alpha=0.85)
        for i, v in enumerate(values):
            ax.text(i, v + 0.3, str(v), ha="center", va="bottom",
                    color="white", fontsize=9, fontweight="bold")
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, fontsize=9, color="#ccccee")
        ax.set_ylim(0, max(values) * 1.3 + 1)
        gen_mark = "ok" if sp.generation_match else "x"
        title = (f"#{idx+1} CY({ba.geometry.h11},{ba.geometry.h21}) "
                 f"-> {sp.gauge_group}\n"
                 f"N_gen={sp.n_generations}({gen_mark})  "
                 f"Score={analysis.final_score:.0f}")
        ax.set_title(title, fontsize=9, color="white", pad=8)

    plt.tight_layout()
    fig.savefig(output_path, dpi=dpi, facecolor=fig.get_facecolor(),
                bbox_inches="tight")
    plt.close(fig)
    return output_path


def plot_sm_comparison(best_analysis, output_path="sm_comparison.png",
                        figsize=(16, 8), dpi=150):
    if best_analysis is None:
        return ""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize, facecolor="#0a0a1a",
                                    gridspec_kw={"width_ratios": [1.2, 1]})
    for ax in (ax1, ax2):
        ax.set_facecolor("#0a0a1a")
        ax.axis("off")

    sp = best_analysis.spectrum
    ba = best_analysis.bundle_analysis

    criteria = [
        ("3 generations", sp.generation_match, sp.n_generations, "= 3"),
        ("Higgs present", sp.higgs_present, sp.n_higgs_candidates, ">= 1"),
        ("Sans exotiques", sp.exotic_free, sp.n_exotics, "= 0"),
        ("c1(V) = 0", ba.c1_vanishes, "-", "oui"),
        ("GUT viable", ba.gut_viable, ba.gauge_group, "SU(5)/SO(10)/E6"),
        ("Anomalie OK", ba.anomaly_cancelled, "-", "oui"),
    ]
    y_start = 0.88
    dy = 0.11
    ax1.text(0.5, 0.98,
             f"Bilan - CY({ba.geometry.h11},{ba.geometry.h21}) -> {sp.gauge_group}",
             fontsize=14, fontweight="bold", color="white",
             ha="center", va="top", transform=ax1.transAxes)

    for i, (name, passed, value, target) in enumerate(criteria):
        y = y_start - i * dy
        color = "#6bcb77" if passed else "#ff6b6b"
        symbol = "+" if passed else "-"
        ax1.text(0.08, y, symbol, fontsize=18, color=color,
                 ha="center", va="center", transform=ax1.transAxes,
                 fontweight="bold")
        ax1.text(0.15, y, name, fontsize=12, color="white",
                 ha="left", va="center", transform=ax1.transAxes)
        ax1.text(0.65, y, str(value), fontsize=12, color=color,
                 ha="center", va="center", transform=ax1.transAxes,
                 fontweight="bold")
        ax1.text(0.88, y, f"(cible: {target})", fontsize=9, color="#888899",
                 ha="center", va="center", transform=ax1.transAxes)

    score = best_analysis.final_score
    sc = "#6bcb77" if score >= 70 else "#ffd93d" if score >= 50 else "#ff6b6b"
    ax1.text(0.5, 0.08, f"Score final : {score:.0f}/100",
             fontsize=18, fontweight="bold", color=sc,
             ha="center", va="center", transform=ax1.transAxes,
             bbox=dict(boxstyle="round,pad=0.5", facecolor="#1a1a2e",
                       edgecolor=sc, alpha=0.8))

    # Panneau droit : representations
    ax2.text(0.5, 0.98, "Decomposition en particules",
             fontsize=13, fontweight="bold", color="white",
             ha="center", va="top", transform=ax2.transAxes)

    if sp.gauge_group == "SU(5)":
        particles = [
            ("10-plet", sp.representations.get("10", 0),
             "q_L, u_R, e_R", "#ff6b6b"),
            ("10bar-plet", sp.representations.get("10bar", 0),
             "anti-particules", "#ff9999"),
            ("5bar-plet", sp.representations.get("5bar", 0),
             "d_R, l_L", "#ffd93d"),
            ("5-plet", sp.representations.get("5", 0),
             "Higgs H_u", "#ffed8a"),
            ("Singlet", sp.representations.get("1", 0),
             "moduli, nu_R", "#8b8b8b"),
        ]
    elif sp.gauge_group == "SO(10)":
        particles = [
            ("16-spineur", sp.representations.get("16", 0),
             "1 generation complete", "#ff6b6b"),
            ("16bar", sp.representations.get("16bar", 0),
             "anti-generation", "#ff9999"),
            ("10-vect", sp.representations.get("10", 0),
             "Higgs", "#ffd93d"),
            ("Singlet", sp.representations.get("1", 0),
             "moduli", "#8b8b8b"),
        ]
    else:
        particles = [
            ("27", sp.representations.get("27", 0),
             "16+10+1 de SO(10)", "#ff6b6b"),
            ("27bar", sp.representations.get("27bar", 0),
             "anti-fondamentale", "#ff9999"),
            ("Singlet", sp.representations.get("1", 0),
             "moduli", "#8b8b8b"),
        ]

    y_p = 0.82
    dy_p = 0.13
    for i, (name, count, desc, color) in enumerate(particles):
        y = y_p - i * dy_p
        bar_w = min(0.8, count * 0.06 + 0.05) if count > 0 else 0.02
        rect = plt.Rectangle((0.05, y - 0.025), bar_w, 0.05,
                              facecolor=color, alpha=0.7, transform=ax2.transAxes)
        ax2.add_patch(rect)
        ax2.text(0.07, y, f"{name}: {count}", fontsize=11,
                 color="white", ha="left", va="center",
                 transform=ax2.transAxes, fontweight="bold")
        ax2.text(bar_w + 0.08, y, desc, fontsize=9, color="#aaaacc",
                 ha="left", va="center", transform=ax2.transAxes)

    plt.tight_layout()
    fig.savefig(output_path, dpi=dpi, facecolor=fig.get_facecolor(),
                bbox_inches="tight")
    plt.close(fig)
    return output_path


def plot_final_ranking(analyses, output_path="final_ranking.png",
                        figsize=(16, 10), dpi=150):
    if not analyses:
        return ""
    n = min(len(analyses), 20)
    top = analyses[:n]
    fig, ax = plt.subplots(figsize=figsize, facecolor="#0a0a1a")
    ax.set_facecolor("#0a0a1a")
    y_positions = range(n)
    bh = 0.35
    bundle_scores = [a.bundle_analysis.sm_score * 0.4 for a in top]
    spectrum_scores = [a.spectrum.sm_compatibility * 0.6 for a in top]

    ax.barh([y - bh/2 for y in y_positions], bundle_scores, height=bh,
            color="#4d96ff", alpha=0.85, label="Fibre (40%)")
    ax.barh([y - bh/2 for y in y_positions], spectrum_scores, height=bh,
            left=bundle_scores, color="#6bcb77", alpha=0.85, label="Spectre (60%)")

    labels = []
    for a in top:
        ba = a.bundle_analysis
        sp = a.spectrum
        gm = "ok" if sp.generation_match else ""
        label = (f"CY({ba.geometry.h11},{ba.geometry.h21}) "
                 f"SU({ba.structure_group_rank})->{sp.gauge_group} "
                 f"[{sp.n_generations}gen{gm}]")
        labels.append(label)

    ax.set_yticks(y_positions)
    ax.set_yticklabels(labels, fontsize=9, color="white", fontfamily="monospace")
    ax.invert_yaxis()

    for i, a in enumerate(top):
        total = a.final_score
        c = "#6bcb77" if total >= 70 else "#ffd93d" if total >= 50 else "#ff6b6b"
        ax.text(total + 1, i, f"{total:.0f}", va="center", ha="left",
                color=c, fontsize=10, fontweight="bold")

    ax.set_xlabel("Score SM combine", color="white", fontsize=12)
    ax.set_xlim(0, 110)
    ax.set_title("Classement final - Pipeline complet CY -> Fibre -> Spectre",
                  fontsize=14, color="white", fontweight="bold", pad=15)
    ax.tick_params(colors="white", labelsize=9)
    for spine in ax.spines.values():
        spine.set_color("#333355")
    ax.legend(loc="lower right", fontsize=10, framealpha=0.3,
              facecolor="#1a1a2e", edgecolor="#333355", labelcolor="white")

    plt.tight_layout()
    fig.savefig(output_path, dpi=dpi, facecolor=fig.get_facecolor(),
                bbox_inches="tight")
    plt.close(fig)
    return output_path


def plot_full_pipeline(best_analysis, output_path="full_pipeline.png",
                        figsize=(18, 7), dpi=150):
    if best_analysis is None:
        return ""
    ba = best_analysis.bundle_analysis
    sp = best_analysis.spectrum

    fig, ax = plt.subplots(figsize=figsize, facecolor="#0a0a1a")
    ax.set_facecolor("#0a0a1a")
    ax.axis("off")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    reps_str = "\n".join(f"{v} x {k}" for k, v in sp.representations.items() if v > 0)

    if sp.gauge_group == "SU(5)":
        sm_p = (f"Quarks   : {sp.n_generations} familles\n"
                f"Leptons  : {sp.n_generations} familles\n"
                f"Higgs    : {sp.n_higgs_candidates} doublet(s)\n"
                f"Exotiques: {sp.n_exotics}")
    elif sp.gauge_group == "SO(10)":
        sm_p = (f"Generations: {sp.n_generations}\n"
                f"(chaque 16 = 1 famille)\n"
                f"Higgs (10) : {sp.n_higgs_candidates}\n"
                f"Exotiques  : {sp.n_exotics}")
    else:
        sm_p = (f"Generations: {sp.n_generations}\n"
                f"(via 27->16+10+1)\n"
                f"Exotiques  : {sp.n_exotics}")

    ck = lambda b: "oui" if b else "non"
    boxes = [
        {"title": "ETAGE 1\nTopologie CY",
         "content": f"h11={ba.geometry.h11}  h21={ba.geometry.h21}\nchi = {ba.geometry.euler}\n|h11-h21| = {ba.geometry.n_gen}",
         "color": "#1a3a6a", "x": 0.02},
        {"title": "ETAGE 2\nFibre vectoriel",
         "content": f"Rang {ba.structure_group_rank} -> SU({ba.structure_group_rank})\nc1=0: {ck(ba.c1_vanishes)}  Anom: {ck(ba.anomaly_cancelled)}\nScore: {ba.sm_score:.0f}/100",
         "color": "#3a1a6a", "x": 0.21},
        {"title": "Groupe\nde jauge",
         "content": f"E8 -> {ba.gauge_group}\nGUT: {ck(ba.gut_viable)}\nBrisure par V",
         "color": "#6a3a1a", "x": 0.40},
        {"title": "ETAGE 3\nSpectre",
         "content": reps_str,
         "color": "#1a6a3a", "x": 0.59},
        {"title": "Particules\ndu SM",
         "content": sm_p,
         "color": "#4a1a4a", "x": 0.78},
    ]

    bw, bh = 0.17, 0.58
    by = 0.24

    for box in boxes:
        rect = plt.Rectangle((box["x"], by), bw, bh,
                              facecolor=box["color"], edgecolor="white",
                              linewidth=1.5, alpha=0.85, transform=ax.transAxes)
        ax.add_patch(rect)
        ax.text(box["x"] + bw/2, by + bh - 0.04, box["title"],
                fontsize=10, fontweight="bold", color="#ffcc66",
                ha="center", va="top", transform=ax.transAxes, linespacing=1.3)
        ax.text(box["x"] + bw/2, by + bh/2 - 0.08, box["content"],
                fontsize=9, color="#ddddee", ha="center", va="center",
                transform=ax.transAxes, linespacing=1.5, fontfamily="monospace")

    arrow_props = dict(arrowstyle="-|>", color="#ffcc66", lw=2.5)
    for i in range(len(boxes) - 1):
        xs = boxes[i]["x"] + bw + 0.003
        xe = boxes[i+1]["x"] - 0.003
        ym = by + bh/2
        ax.annotate("", xy=(xe, ym), xytext=(xs, ym),
                     arrowprops=arrow_props, transform=ax.transAxes)

    score = best_analysis.final_score
    sc = "#6bcb77" if score >= 70 else "#ffd93d" if score >= 50 else "#ff6b6b"
    ax.text(0.5, 0.08, f"Score final combine : {score:.0f}/100",
            fontsize=16, fontweight="bold", color=sc,
            ha="center", va="center", transform=ax.transAxes,
            bbox=dict(boxstyle="round,pad=0.6", facecolor="#0a0a1a",
                      edgecolor=sc, linewidth=2, alpha=0.9))

    ax.set_title("Pipeline complet : de la geometrie aux particules elementaires",
                 fontsize=16, color="white", fontweight="bold", y=0.96)

    fig.savefig(output_path, dpi=dpi, facecolor=fig.get_facecolor(),
                bbox_inches="tight")
    plt.close(fig)
    return output_path
