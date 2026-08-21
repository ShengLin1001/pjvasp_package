#!/usr/bin/env python3
"""Generate post-processing documentation figures from real data.

Reads data files from /public3/home/scg6928/mywork/temp/post-docs/data/
Saves PNG figures to /public3/home/scg6928/mywork/temp/post-docs/figures/

Uses mymetal plotting style: general_set_all_rcParams + my_plot + general_modify_legend.
Run with: MPLBACKEND=Agg <codex_python> gen_post_docs_figures.py
"""
import os
import re
import json
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

# --- mymetal style imports (use general_set_all_rcParams as the entry point) ---
from mymetal.universal.plot.general import general_set_all_rcParams, general_modify_legend
from mymetal.universal.plot.plot import my_plot

# ── Paths ──────────────────────────────────────────────────────────────────
DATA_DIR = "/public3/home/scg6928/mywork/temp/post-docs/data"
FIG_DIR  = "/public3/home/scg6928/mywork/temp/post-docs/figures"
os.makedirs(FIG_DIR, exist_ok=True)

# ── Set global style once ──────────────────────────────────────────────────
general_set_all_rcParams()


def _save(fig, name):
    path = os.path.join(FIG_DIR, name)
    fig.savefig(path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"  saved: {path}  ({os.path.getsize(path)} bytes)")


# ════════════════════════════════════════════════════════════════════════════
#  Data loaders  (parsers for the various text/JSON files)
# ════════════════════════════════════════════════════════════════════════════

def load_hoec():
    with open(os.path.join(DATA_DIR, "hoec_au_hex.json")) as f:
        return json.load(f)

def load_cij_vasp():
    """Parse the Cij block from cij_au_vasp.txt — returns dict of {label: value}."""
    txt = open(os.path.join(DATA_DIR, "cij_au_vasp.txt")).read()
    # first data row after 'C11 … C44' header
    m = re.search(r'C11\s+C12\s+C13\s+C33\s+C44\s*\n\s*([\d.\s]+)', txt)
    vals = [float(x) for x in m.group(1).split()]
    return {"C11": vals[0], "C12": vals[1], "C13": vals[2], "C33": vals[3], "C44": vals[4]}

def load_cij_lammps():
    txt = open(os.path.join(DATA_DIR, "cij_fcc_lammps.txt")).read()
    m = re.search(r'C11\s+C12\s+C13\s+C33\s+C44\s*\n\s*([\d.\s]+)', txt)
    vals = [float(x) for x in m.group(1).split()]
    return {"C11": vals[0], "C12": vals[1], "C13": vals[2], "C33": vals[3], "C44": vals[4]}

def load_stretch_vasp():
    """Return (factors, energies) arrays from the VASP stretch file."""
    txt = open(os.path.join(DATA_DIR, "stretch_au_vasp.txt")).read()
    factors, energies = [], []
    for line in txt.splitlines():
        parts = line.split()
        if len(parts) == 3:
            try:
                f = float(parts[0])
                e = float(parts[1])
                # data rows have factor near 1.0
                if 0.99 < f < 1.01:
                    factors.append(f)
                    energies.append(e)
            except ValueError:
                pass
    return np.array(factors), np.array(energies)

def load_stretch_lammps():
    txt = open(os.path.join(DATA_DIR, "stretch_fcc_lammps.txt")).read()
    factors, energies = [], []
    for line in txt.splitlines():
        parts = line.split()
        if len(parts) == 3:
            try:
                f = float(parts[0])
                e = float(parts[1])
                if 0.99 < f < 1.01:
                    factors.append(f)
                    energies.append(e)
            except ValueError:
                pass
    return np.array(factors), np.array(energies)

def load_gsfe():
    """Return (displacement_norm, gamma) from GSFE file — displacement is slip column (Ang)."""
    txt = open(os.path.join(DATA_DIR, "gsfe_fcc_111_lammps.txt")).read()
    # columns: jobn  dE  gamma  da31/a11  da32/a22  slip  da33
    gamma = []
    slip = []
    for line in txt.splitlines():
        parts = line.split()
        if len(parts) >= 6:
            try:
                g = float(parts[2])
                s = float(parts[5])
                gamma.append(g)
                slip.append(s)
            except ValueError:
                pass
    # also extract gamma_usf (local max)
    m = re.search(r'local max \(mJ/m\^2\):\s*([\d.]+)', txt)
    gamma_usf = float(m.group(1)) if m else max(gamma)
    return np.array(slip), np.array(gamma), gamma_usf

def load_encut():
    txt = open(os.path.join(DATA_DIR, "convergence_encuts_au_fcc.txt")).read()
    encuts, energies = [], []
    for line in txt.splitlines():
        parts = line.split()
        if len(parts) == 2:
            try:
                e = float(parts[0])      # encut (eV)
                et = float(parts[1])    # Etot
                if 200 < e < 700:
                    encuts.append(e)
                    energies.append(et)
            except ValueError:
                pass
    return np.array(encuts), np.array(energies)

def load_kpoints():
    txt = open(os.path.join(DATA_DIR, "convergence_kpoints_au_fcc.txt")).read()
    kgrids, energies = [], []
    for line in txt.splitlines():
        parts = line.split()
        if len(parts) == 2:
            kstr = parts[0]
            try:
                et = float(parts[1])
                # parse '10-10-10' → 10
                nums = re.findall(r'(\d+)-(\d+)-(\d+)', kstr)
                if nums:
                    n = int(nums[0][0])
                    kgrids.append(n)
                    energies.append(et)
            except (ValueError, IndexError):
                pass
    # sort by kgrid
    order = np.argsort(kgrids)
    return np.array(kgrids)[order], np.array(energies)[order]

def load_relax():
    txt = open(os.path.join(DATA_DIR, "relax_convergence_00.txt")).read()
    frames, energies, forces = [], [], []
    for line in txt.splitlines():
        parts = line.split()
        if len(parts) == 3:
            try:
                fr = int(float(parts[0]))
                en = float(parts[1])
                fo = float(parts[2])
                frames.append(fr)
                energies.append(en)
                forces.append(fo)
            except ValueError:
                pass
    return np.array(frames), np.array(energies), np.array(forces)


# ════════════════════════════════════════════════════════════════════════════
#  Figure 1: post_hoec_energy_au.png  —  HOEC bar chart (2nd/3rd/4th order)
# ════════════════════════════════════════════════════════════════════════════

def fig1_hoec_bar():
    print("[1/8] HOEC bar chart ...")
    data = load_hoec()
    c2 = data["constants"]["second_order"]
    c3 = data["constants"]["third_order"]
    c4 = data["constants"]["fourth_order"]

    # Grouped bar: 2nd (5 constants), 3rd (10), 4th (19)
    labels2 = ["C11", "C12", "C13", "C33", "C44"]
    labels3 = ["C111", "C112", "C113", "C123", "C133", "C144", "C155", "C222", "C333", "C344"]
    labels4 = ["C1111", "C1112", "C1113", "C1122", "C1123", "C1133", "C1144", "C1155",
               "C1166", "C1223", "C1233", "C1244", "C1255", "C1333", "C1344", "C1355",
               "C3333", "C3344", "C4444"]

    vals2 = [float(c2[l]) for l in labels2]
    vals3 = [float(c3[l]) for l in labels3]
    vals4 = [float(c4[l]) for l in labels4]

    fig, axes = my_plot(fig_subp=[3, 1], fig_sharex=False)
    axes = axes.flatten() if hasattr(axes, 'flatten') else axes

    colors = ["#2166ac", "#b2182b", "#4daf4a"]
    titles = ["2nd order (SOEC)", "3rd order (TOEC)", "4th order (FOEC)"]
    all_data = [(labels2, vals2), (labels3, vals3), (labels4, vals4)]

    for i, (ax, (labs, vals), title, color) in enumerate(zip(axes, all_data, titles, colors)):
        x = np.arange(len(labs))
        ax.bar(x, vals, color=color, edgecolor='black', linewidth=1.5, width=0.7)
        ax.set_xticks(x)
        ax.set_xticklabels(labs, rotation=45, ha='right', fontsize=16)
        ax.set_ylabel("Elastic constant (GPa)")
        ax.set_title(title, fontsize=22)
        # log scale for magnitude comparison
        ax.set_yscale('symlog', linthresh=10)
        ax.axhline(0, color='gray', linewidth=1, linestyle='--')

    fig.savefig(os.path.join(FIG_DIR, "post_hoec_energy_au.png"),
                dpi=300, bbox_inches='tight')
    plt.close(fig)
    sz = os.path.getsize(os.path.join(FIG_DIR, "post_hoec_energy_au.png"))
    print(f"  saved: post_hoec_energy_au.png ({sz} bytes)")


# ════════════════════════════════════════════════════════════════════════════
#  Figure 2: post_cij_comparison.png  —  VASP-HCP vs LAMMPS-FCC Cij bars
# ════════════════════════════════════════════════════════════════════════════

def fig2_cij_comparison():
    print("[2/8] Cij comparison bar chart ...")
    cij_v = load_cij_vasp()
    cij_l = load_cij_lammps()

    labels = ["C11", "C12", "C13", "C33", "C44"]
    vals_v = [float(cij_v[l]) for l in labels]
    vals_l = [float(cij_l[l]) for l in labels]

    x = np.arange(len(labels))
    w = 0.35

    fig, ax = my_plot()
    ax = ax[0] if hasattr(ax, '__getitem__') and not hasattr(ax, 'plot') else ax

    b1 = ax.bar(x - w/2, vals_v, w, label="VASP (Au HCP)", color="#2166ac",
                edgecolor='black', linewidth=1.5)
    b2 = ax.bar(x + w/2, vals_l, w, label="LAMMPS (FCC)", color="#d6604d",
                edgecolor='black', linewidth=1.5)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=24)
    ax.set_ylabel("Elastic constant (GPa)")
    ax.set_title("Elastic constants: VASP (Au HCP) vs LAMMPS (FCC)", fontsize=22)
    general_modify_legend(ax.legend())

    fig.savefig(os.path.join(FIG_DIR, "post_cij_comparison.png"),
                dpi=300, bbox_inches='tight')
    plt.close(fig)
    sz = os.path.getsize(os.path.join(FIG_DIR, "post_cij_comparison.png"))
    print(f"  saved: post_cij_comparison.png ({sz} bytes)")


# ════════════════════════════════════════════════════════════════════════════
#  Figure 3: post_stretch_comparison.png  —  two-panel stretch curves w/ fit
# ════════════════════════════════════════════════════════════════════════════

def fig3_stretch_comparison():
    print("[3/8] Stretch comparison (two-panel) ...")
    fv, ev = load_stretch_vasp()
    fl, el = load_stretch_lammps()

    # per-atom energies
    ev_pa = ev / 2.0   # VASP natoms=2
    el_pa = el / 4.0   # LAMMPS natoms=4

    # quadratic fit: E = a*(f-1)^2 + b*(f-1) + c
    dv = (fv - 1.0)
    dl = (fl - 1.0)
    # shift energy so minimum is 0 for readability
    ev_rel = ev_pa - ev_pa.min()
    el_rel = el_pa - el_pa.min()

    cv = np.polyfit(dv, ev_rel, 2)
    cl = np.polyfit(dl, el_rel, 2)

    fig, axes = my_plot(fig_subp=[1, 2], fig_sharex=False)
    ax1, ax2 = axes[0], axes[1]

    # VASP panel
    ax1.plot(fv, ev_rel * 1000, marker='o', color="#2166ac", linewidth=3,
             markersize=10, markerfacecolor='white', label="VASP data")
    xx = np.linspace(fv.min(), fv.max(), 200)
    ax1.plot(xx, np.polyval(cv, xx - 1.0) * 1000, '--', color="#2166ac",
             linewidth=2, label="Quadratic fit")
    ax1.set_xlabel("Stretch factor")
    ax1.set_ylabel(r"$\Delta E$ (meV/atom)")
    ax1.set_title("VASP — Au HCP (xy stretch)", fontsize=20)
    general_modify_legend(ax1.legend())

    # LAMMPS panel
    ax2.plot(fl, el_rel * 1000, marker='s', color="#d6604d", linewidth=3,
             markersize=8, markerfacecolor='white', label="LAMMPS data")
    xx = np.linspace(fl.min(), fl.max(), 200)
    ax2.plot(xx, np.polyval(cl, xx - 1.0) * 1000, '--', color="#d6604d",
             linewidth=2, label="Quadratic fit")
    ax2.set_xlabel("Stretch factor")
    ax2.set_ylabel(r"$\Delta E$ (meV/atom)")
    ax2.set_title("LAMMPS — FCC (xyz stretch)", fontsize=20)
    general_modify_legend(ax2.legend())

    fig.savefig(os.path.join(FIG_DIR, "post_stretch_comparison.png"),
                dpi=300, bbox_inches='tight')
    plt.close(fig)
    sz = os.path.getsize(os.path.join(FIG_DIR, "post_stretch_comparison.png"))
    print(f"  saved: post_stretch_comparison.png ({sz} bytes)")


# ════════════════════════════════════════════════════════════════════════════
#  Figure 4: post_gsfe_fcc_111.png  —  GSFE curve with SF/USF markers
# ════════════════════════════════════════════════════════════════════════════

def fig4_gsfe():
    print("[4/8] GSFE FCC-111 curve ...")
    slip, gamma, gamma_usf = load_gsfe()

    fig, ax = my_plot()
    ax = ax[0] if hasattr(ax, '__getitem__') and not hasattr(ax, 'plot') else ax

    ax.plot(slip, gamma, marker='o', color="#2166ac", linewidth=3, markersize=12,
            markerfacecolor='white', label="GSFE curve")

    # Mark USF (unstable stacking fault) — the local maximum
    idx_usf = int(np.argmax(gamma))
    ax.annotate(f"$\\gamma_{{usf}}$ = {round(float(gamma_usf), 1)} mJ/m$^2$",
                xy=(slip[idx_usf], gamma[idx_usf]),
                xytext=(slip[idx_usf] + 0.3, gamma[idx_usf] + 5),
                fontsize=20, arrowprops=dict(arrowstyle='->', color='red', lw=2),
                color='red')
    ax.plot(slip[idx_usf], gamma[idx_usf], 'r*', markersize=20, label="USF (local max)")

    # Mark SF (stable stacking fault) — local minimum after the peak, at ~2/3 displacement
    # The stable SF is at the minimum after USF
    after_peak = gamma[idx_usf + 1:]
    if len(after_peak) > 0:
        idx_sf_local = int(np.argmin(after_peak))
        idx_sf = idx_usf + 1 + idx_sf_local
        gamma_sf = float(gamma[idx_sf])
        ax.plot(slip[idx_sf], gamma[idx_sf], 'g^', markersize=18, label=f"SF = {round(gamma_sf, 1)} mJ/m$^2$")

    ax.set_xlabel("Displacement along $\\langle 11\\bar{2} \\rangle$ (Å)")
    ax.set_ylabel(r"$\gamma$ (mJ/m$^2$)")
    ax.set_title("Generalized Stacking Fault Energy — FCC (111)", fontsize=22)
    general_modify_legend(ax.legend())

    fig.savefig(os.path.join(FIG_DIR, "post_gsfe_fcc_111.png"),
                dpi=300, bbox_inches='tight')
    plt.close(fig)
    sz = os.path.getsize(os.path.join(FIG_DIR, "post_gsfe_fcc_111.png"))
    print(f"  saved: post_gsfe_fcc_111.png ({sz} bytes)")


# ════════════════════════════════════════════════════════════════════════════
#  Figure 5: post_convergence_encuts.png  —  ENCUT convergence
# ════════════════════════════════════════════════════════════════════════════

def fig5_encut_conv():
    print("[5/8] ENCUT convergence ...")
    encuts, energies = load_encut()

    # converged reference = highest encut energy
    e_conv = float(energies[-1])
    delta = (energies - e_conv) * 1000  # meV

    fig, ax = my_plot()
    ax = ax[0] if hasattr(ax, '__getitem__') and not hasattr(ax, 'plot') else ax

    ax.plot(encuts, delta, marker='o', color="#2166ac", linewidth=3, markersize=14,
            markerfacecolor='white', label=r"$\Delta E$ from converged")
    ax.axhline(0, color='gray', linestyle='--', linewidth=1.5)
    # annotate converged point
    ax.annotate(f"Converged\n({int(encuts[-1])} eV)",
                xy=(encuts[-1], delta[-1]),
                xytext=(encuts[-1] - 100, delta[-1] + 5),
                fontsize=18, arrowprops=dict(arrowstyle='->', color='red', lw=2),
                color='red')

    ax.set_xlabel("ENCUT (eV)")
    ax.set_ylabel(r"$\Delta E$ from converged (meV)")
    ax.set_title("Energy convergence vs ENCUT (Au FCC)", fontsize=22)
    general_modify_legend(ax.legend())

    fig.savefig(os.path.join(FIG_DIR, "post_convergence_encuts.png"),
                dpi=300, bbox_inches='tight')
    plt.close(fig)
    sz = os.path.getsize(os.path.join(FIG_DIR, "post_convergence_encuts.png"))
    print(f"  saved: post_convergence_encuts.png ({sz} bytes)")


# ════════════════════════════════════════════════════════════════════════════
#  Figure 6: post_convergence_kpoints.png  —  KPOINTS convergence
# ════════════════════════════════════════════════════════════════════════════

def fig6_kpoints_conv():
    print("[6/8] KPOINTS convergence ...")
    kgrids, energies = load_kpoints()

    e_conv = float(energies[-1])
    delta = (energies - e_conv) * 1000  # meV

    fig, ax = my_plot()
    ax = ax[0] if hasattr(ax, '__getitem__') and not hasattr(ax, 'plot') else ax

    ax.plot(kgrids, delta, marker='o', color="#b2182b", linewidth=3, markersize=12,
            markerfacecolor='white', label=r"$\Delta E$ from converged")
    ax.axhline(0, color='gray', linestyle='--', linewidth=1.5)

    ax.set_xlabel("K-grid ($N \\times N \\times N$)")
    ax.set_ylabel(r"$\Delta E$ from converged (meV)")
    ax.set_title("Energy convergence vs K-mesh (Au FCC)", fontsize=22)
    general_modify_legend(ax.legend())

    fig.savefig(os.path.join(FIG_DIR, "post_convergence_kpoints.png"),
                dpi=300, bbox_inches='tight')
    plt.close(fig)
    sz = os.path.getsize(os.path.join(FIG_DIR, "post_convergence_kpoints.png"))
    print(f"  saved: post_convergence_kpoints.png ({sz} bytes)")


# ════════════════════════════════════════════════════════════════════════════
#  Figure 7: post_relax_convergence.png  —  two-panel energy/force vs ionic step
# ════════════════════════════════════════════════════════════════════════════

def fig7_relax():
    print("[7/8] Relaxation convergence ...")
    frames, energies, forces = load_relax()

    fig, axes = my_plot(fig_subp=[1, 2], fig_sharex=False)
    ax1, ax2 = axes[0], axes[1]

    # Energy panel — shift to min for readability
    e_rel = energies - energies.min()
    ax1.plot(frames, e_rel * 1000, marker='o', color="#2166ac", linewidth=3,
             markersize=14, markerfacecolor='white')
    ax1.set_xlabel("Ionic step")
    ax1.set_ylabel(r"$\Delta E$ from min (meV)")
    ax1.set_title("Energy convergence", fontsize=20)

    # Force panel — log scale (forces are 0 in this data, handle gracefully)
    if forces.max() > 0:
        ax2.semilogy(frames, forces, marker='s', color="#d6604d", linewidth=3,
                     markersize=14, markerfacecolor='white')
    else:
        ax2.plot(frames, forces, marker='s', color="#d6604d", linewidth=3,
                 markersize=14, markerfacecolor='white')
    ax2.set_xlabel("Ionic step")
    ax2.set_ylabel("Max force (eV/Å)")
    ax2.set_title("Force convergence", fontsize=20)

    fig.savefig(os.path.join(FIG_DIR, "post_relax_convergence.png"),
                dpi=300, bbox_inches='tight')
    plt.close(fig)
    sz = os.path.getsize(os.path.join(FIG_DIR, "post_relax_convergence.png"))
    print(f"  saved: post_relax_convergence.png ({sz} bytes)")


# ════════════════════════════════════════════════════════════════════════════
#  Figure 8: post_hoec_mode_fits.png  —  per-mode energy-strain fits from HOEC
# ════════════════════════════════════════════════════════════════════════════

def fig8_hoec_mode_fits():
    print("[8/8] HOEC per-mode energy-strain fits ...")
    data = load_hoec()

    # We have the solved constants.  Reconstruct per-mode energy curves using
    # E = 1/2*C2*η² + 1/6*C3*η³ + 1/24*C4*η⁴ for a representative uniaxial
    # mode.  Use C11 family (C11, C111, C1111) as a proxy "M01" mode.
    c2 = float(data["constants"]["second_order"]["C11"])
    c3 = float(data["constants"]["third_order"]["C111"])
    c4 = float(data["constants"]["fourth_order"]["C1111"])

    # strain range matching fitmax
    fitmax = float(data.get("fitmax", 0.12))
    eta = np.linspace(-fitmax, fitmax, 200)

    e2 = 0.5 * c2 * eta**2
    e3 = (1.0 / 6.0) * c3 * eta**3
    e4 = (1.0 / 24.0) * c4 * eta**4
    e_tot = e2 + e3 + e4

    fig, ax = my_plot()
    ax = ax[0] if hasattr(ax, '__getitem__') and not hasattr(ax, 'plot') else ax

    ax.plot(eta, e_tot, '-', color="black", linewidth=3, label="Total (2+3+4)")
    ax.plot(eta, e2, '--', color="#2166ac", linewidth=2.5, label="2nd order")
    ax.plot(eta, e2 + e3, '--', color="#d6604d", linewidth=2.5, label="2nd+3rd")
    ax.plot(eta, e2 + e3 + e4, ':', color="#4daf4a", linewidth=2.5,
            label="2nd+3rd+4th")

    ax.set_xlabel("Strain η")
    ax.set_ylabel("Energy density (GPa)")
    ax.set_title("HOEC mode fit — C11 family (Au HCP)", fontsize=22)
    general_modify_legend(ax.legend())

    fig.savefig(os.path.join(FIG_DIR, "post_hoec_mode_fits.png"),
                dpi=300, bbox_inches='tight')
    plt.close(fig)
    sz = os.path.getsize(os.path.join(FIG_DIR, "post_hoec_mode_fits.png"))
    print(f"  saved: post_hoec_mode_fits.png ({sz} bytes)")


# ════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    print(f"Data dir: {DATA_DIR}")
    print(f"Figure dir: {FIG_DIR}")
    print()
    fig1_hoec_bar()
    fig2_cij_comparison()
    fig3_stretch_comparison()
    fig4_gsfe()
    fig5_encut_conv()
    fig6_kpoints_conv()
    fig7_relax()
    fig8_hoec_mode_fits()
    print("\n=== All figures generated ===")
    # list all PNGs
    for f in sorted(os.listdir(FIG_DIR)):
        if f.endswith('.png'):
            sz = os.path.getsize(os.path.join(FIG_DIR, f))
            status = "OK" if sz > 5000 else "WARN (<5KB)"
            print(f"  {f:40s}  {sz:>8d} bytes  [{status}]")
