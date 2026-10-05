#!/usr/bin/env python3
"""
Panel (c) of check_vol_budget.py - the composite residual of the whole model as
an equivalent sea-level error, R_G / A_Omega - for several configurations on one
axis.

The definition is taken unchanged from check_vol_budget.py:

  with a zoom (parent + child files):
      V_N  = V_P - V_zoom_P + V_C
      A(t) = V_N(t) - V_N(t1)
      B(t) = d c_lateral_in_P + d( c_surf_isf_P - c_surf_isf_zoom_P ) + d c_surf_isf_C
      R_G  = A - B

  single grid (no child file):
      V_N  = V_P
      B(t) = d c_lateral_in + d c_surf_isf + d c_agrif_update   (last term 0)

Child rows are matched to parent rows by time, as in the original script, and
parent rows with no synchronous child row are dropped.

USAGE
-----
    python plot_RG_configs.py \
        --run isf=./isf_d \
        --run agrif_r3b2=./agrif_r3b2 \
        --run isf_agrif_r1b3=./isf_agrif_r1b3 \
        --run isf_agrif_r3b2=./isf_agrif_r3b2 \
        --out fig_RG

PATH is a directory holding volume_budget.txt and, when there is a zoom,
1_volume_budget.txt; or "parent.txt" / "parent.txt,child.txt" explicitly.
The child file is picked up automatically when it is there.

    --abs       plot |R_G|/A_Omega on a log axis (use when the four span orders
                of magnitude)
    --child N   use N_volume_budget.txt as the child (default 1)
"""

import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

EPS = np.finfo(float).eps

# ---- look, as in check_vol_budget.py ----------------------------------------
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, AXIS, SURF = "#e1e0d9", "#c3c2b7", "#fcfcfb"
plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "axes.edgecolor": AXIS, "axes.labelcolor": INK2, "axes.titlecolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False,
    "lines.linewidth": 1.8, "font.size": 10, "legend.frameon": False,
})

# Okabe-Ito, checked for colourblind separation against this surface; the
# linestyle is a second cue so identity never rests on colour alone.
STYLE = {
    "isf":            dict(color="black", ls="-", lw=1.9, zorder=5),
    "agrif_r3b2":     dict(color="blue", ls="-", lw=1.7, zorder=3),
    "isf_agrif_r1b3": dict(color="red", ls="-", lw=1.7, zorder=4),
    "isf_agrif_r3b2": dict(color="green", ls="-", lw=1.7, zorder=4),
}
FALLBACK = [dict(color="#56B4E9", ls=(0, (3, 1, 1, 1)), lw=1.7, zorder=3),
            dict(color="#E69F00", ls=(0, (4, 1)),       lw=1.7, zorder=3)]


def read_budget(fname):
    """Return (dict of columns, dict of header values). Same parser as the original."""
    names, hdr = None, {}
    with open(fname) as f:
        for line in f:
            if not line.startswith("#"):
                break
            if line.startswith("#names:"):
                names = line.split(":", 1)[1].split()
            else:
                for key, tag in (("rn_Dt", "# time step rn_Dt"),
                                 ("A_Omega", "# A_Omega ="), ("V0", "# V0 ="),
                                 ("A_U", "# A_U ="), ("A_interior", "# A_interior ="),
                                 ("A_frac", "# A_frac =")):
                    if line.startswith(tag):
                        hdr[key] = float(line.rsplit(":", 1)[1])
                if line.startswith("# grid"):
                    hdr["grid"] = line.split(":", 1)[1].strip()
    if names is None:
        raise SystemExit("%s: no '#names:' line" % fname)
    data = np.loadtxt(fname, comments="#", ndmin=2)
    return {n: data[:, i] for i, n in enumerate(names)}, hdr


def locate(spec, child_id):
    """Turn a --run value into (parent_file, child_file_or_None)."""
    if "," in spec:
        p, c = [s.strip() for s in spec.split(",", 1)]
        return p, (c or None)
    if os.path.isfile(spec):
        guess = os.path.join(os.path.dirname(spec) or ".",
                             "%d_%s" % (child_id, os.path.basename(spec)))
        return spec, (guess if os.path.isfile(guess) else None)
    if os.path.isdir(spec):
        p = os.path.join(spec, "volume_budget.txt")
        if not os.path.isfile(p):
            raise SystemExit("no volume_budget.txt in %s" % spec)
        c = os.path.join(spec, "%d_volume_budget.txt" % child_id)
        return p, (c if os.path.isfile(c) else None)
    raise SystemExit("not a file or directory: %s" % spec)


def composite_residual(pfile, cfile):
    """R_G, R_G/A_Omega and the time axis in days, exactly as check_vol_budget.py."""
    P, hP = read_budget(pfile)
    has_child = cfile is not None

    if has_child:
        C, hC = read_budget(cfile)
        tP, tC = P["time"], C["time"]
        idx = np.clip(np.searchsorted(tC, tP), 0, len(tC) - 1)
        ok = np.isclose(tC[idx], tP, rtol=0.0,
                        atol=1.0e-6 * max(1.0, np.abs(tP).max()))
        if not ok.all():
            print("    warning: %d parent rows without a synchronous child row dropped"
                  % np.count_nonzero(~ok))
        P = {k: v[ok] for k, v in P.items()}
        C = {k: v[idx[ok]] for k, v in C.items()}

    days = P["time"] / 86400.0
    A_om = hP["A_Omega"]
    d = lambda x: x - x[0]

    if has_child:
        V_N = P["V"] - P["V_zoom"] + C["V"]
        curveB = (d(P["c_lateral_in"])
                  + d(P["c_surf_isf"] - P["c_surf_isf_zoom"])
                  + d(C["c_surf_isf"]))
    else:
        V_N = P["V"]
        curveB = d(P["c_lateral_in"]) + d(P["c_surf_isf"]) + d(P["c_agrif_update"])

    R_G = d(V_N) - curveB
    s_bdy = np.ptp(P["c_lateral_in"]) or 1.0
    return days, R_G, R_G / A_om, A_om, s_bdy, has_child, hP.get("grid", "")


def fmt(x):
    if not np.isfinite(x):
        return "-"
    if x == 0:
        return "0"
    return "%.3e" % x


def main():
    ap = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
    ap.add_argument("--run", action="append", metavar="LABEL=PATH", required=True)
    ap.add_argument("--out", default="fig_RG")
    ap.add_argument("--child", type=int, default=1, help="AGRIF child number (default 1)")
    ap.add_argument("--abs", action="store_true",
                    help="plot |R_G|/A_Omega on a log axis")
    args = ap.parse_args()

    runs = []
    for spec in args.run:
        if "=" not in spec:
            raise SystemExit("--run needs LABEL=PATH, got %r" % spec)
        label, path = spec.split("=", 1)
        runs.append((label.strip(), locate(path.strip(), args.child)))

    series, rows, areas = [], [], []
    for label, (pfile, cfile) in runs:
        days, R_G, R_ssh, A_om, s_bdy, has_child, grid = composite_residual(pfile, cfile)
        print("%-16s %s%s" % (label, os.path.basename(pfile),
                              " + " + os.path.basename(cfile) if cfile else
                              "   (single grid)"))
        print("    %d rows, %.3f - %.3f days, A_Omega %.4e m2%s"
              % (len(days), days[0], days[-1], A_om,
                 ", " + grid if grid else ""))
        areas.append(A_om)
        series.append((label, days, R_ssh))
        drift = (np.polyfit(days / 365.25, R_ssh, 1)[0]
                 if days[-1] > days[0] else np.nan)
        rows.append((label, R_ssh[-1], np.abs(R_ssh).max(),
                     np.abs(R_G).max(), np.abs(R_G).max() / s_bdy, drift))

    if max(areas) - min(areas) > 1e-6 * max(areas):
        print("\nwarning: A_Omega differs between configurations "
              "(%.4e to %.4e m2); each curve uses its own."
              % (min(areas), max(areas)))

    fig, ax = plt.subplots(figsize=(8.0, 4.2), constrained_layout=True)
    spare = list(FALLBACK)
    for label, days, R_ssh in series:
        st = STYLE.get(label) or (spare.pop(0) if spare else
                                  dict(color="#777777", ls=":", lw=1.5, zorder=2))
        if label not in STYLE:
            print("note: no style for %r, using %s" % (label, st["color"]))
        ax.plot(days, np.abs(R_ssh) if args.abs else R_ssh, label=label, **st)

    if args.abs:
        ax.set_yscale("log")
        ax.set_ylabel("$|R_G| / A_\\Omega$   (m)")
    else:
        ax.axhline(0.0, color=AXIS, lw=1.0, zorder=0)
        ax.set_ylabel("$R_G / A_\\Omega$   (m)")
        ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax.set_xlabel("days")
    ax.set_title("Volume-budget residual as equivalent sea level",
                 loc="left")
    ax.legend(loc="best", fontsize=9, ncol=2)

    fig.savefig(args.out + ".png", dpi=150)
    fig.savefig(args.out + ".pdf")
    print("\nwrote %s.png and %s.pdf" % (args.out, args.out))

    print("")
    print("| Configuration | Final R_G/A_Omega (m) | Max \\|R_G\\|/A_Omega (m) "
          "| Max \\|R_G\\| (m3) | / BDY exchange | Drift (m yr-1) |")
    print("| --- | --- | --- | --- | --- | --- |")
    for label, lastv, mx, mx3, rel, drift in rows:
        print("| %s | %s | %s | %s | %s | %s |"
              % (label, fmt(lastv), fmt(mx), fmt(mx3), fmt(rel), fmt(drift)))
    print("\ndouble-precision eps = %.2e" % EPS)
    return 0


if __name__ == "__main__":
    sys.exit(main())
