"""a matplotlib dashboard for the anafield benchmark inputs.

sliders for the grid size, the variogram range, the nonlinearity a, the
observation counts, the window size, the noise and the seed.  the panels show
the truth field with the observation layout drawn on it, and for each
observation type the nonlinear observation against its linear counterpart so
you can see how much a, the path length and the window size bend things.
the write button calls anafield_setup.setup with the slider values.

    python anafield_dashboard.py
"""
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from matplotlib.gridspec import GridSpec

import anafield_setup as af


class Dashboard:

    def __init__(self, out_d="anafield_dash_template"):
        self.out_d = out_d
        self.fig = plt.figure(figsize=(15, 9.5))
        self.fig.canvas.manager.set_window_title("anafield benchmark inputs")
        gs = GridSpec(2, 3, left=0.04, right=0.70, bottom=0.07, top=0.95, wspace=0.3, hspace=0.3)
        self.ax_field = self.fig.add_subplot(gs[:, :2])
        self.ax_point = self.fig.add_subplot(gs[0, 2])
        self.ax_path = self.fig.add_subplot(gs[1, 2])
        self.sliders = {}
        x0, w, h, dy = 0.78, 0.18, 0.018, 0.034
        y = 0.92
        self.fig.text(x0, y + 0.03, "inputs", fontsize=11, weight="bold")

        def add(key, label, lo, hi, val, step=None, fmt=None):
            nonlocal y
            ax = self.fig.add_axes([x0, y, w, h])
            s = Slider(ax, label, lo, hi, valinit=val, valstep=step, valfmt=fmt)
            s.label.set_fontsize(8.5)
            s.valtext.set_fontsize(8.5)
            s.on_changed(self.schedule)
            self.sliders[key] = s
            y -= dy

        self.fig.text(x0, y + 0.005, "field", fontsize=8.5, style="italic")
        y -= 0.03
        add("nrow", "nrow", 5, 200, 20, 1, "%d")
        add("ncol", "ncol", 5, 200, 20, 1, "%d")
        add("range", "range (cells)", 1, 50, 8, 1, "%d")
        add("sill", "log10 variance", 0.01, 1.0, 0.25, 0.01, "%.2f")
        add("seed", "seed", 1, 500, 12345 % 500, 1, "%d")
        y -= 0.015
        self.fig.text(x0, y + 0.005, "observations", fontsize=8.5, style="italic")
        y -= 0.03
        add("a", "a (point nonlinearity)", 0.0, 4.0, 2.0, 0.1, "%.1f")
        add("n_points", "point obs", 0, 200, 20, 1, "%d")
        add("n_paths", "path obs", 0, 100, 10, 1, "%d")
        add("n_windows", "block obs", 0, 100, 5, 1, "%d")
        add("window", "window (cells)", 1, 40, 4, 1, "%d")
        add("noise", "noise std (log10)", 0.005, 0.5, 0.05, 0.005, "%.3f")
        add("num_reals", "prior reals", 10, 500, 100, 10, "%d")

        y -= 0.02
        ax = self.fig.add_axes([x0, y - 0.01, 0.08, 0.035])
        self.btn = Button(ax, "write setup")
        self.btn.label.set_fontsize(8.5)
        self.btn.on_clicked(self.write)
        self.status = self.fig.text(x0 + 0.09, y + 0.005, "", fontsize=7.5, va="center")
        self.summary = self.fig.text(x0, y - 0.045, "", fontsize=8, va="top", family="monospace")

        defs = [
            "nrow, ncol   grid size; parameters = nrow x ncol",
            "range        variogram e-folding length in cells;",
            "             the practical range is about 3x",
            "log10 var    sill; 0.25 = log10 std of 0.5",
            "a            point obs = sinh(a x)/a;",
            "             0 = linear, 3 = strongly convex",
            "path obs     log10 harmonic mean of 10^x along",
            "             a random line: a soft minimum",
            "block obs    log10 arithmetic mean of 10^x over",
            "             a window: a soft maximum",
            "window       side of the block, in cells",
            "noise        std added to every obs; weight 1/noise",
            "prior reals  ensemble size written to prior_pe.jcb",
            "write setup  calls anafield_setup.setup() with these",
        ]
        self.fig.text(x0, 0.27, "what the controls mean", fontsize=9, weight="bold", va="top")
        self.fig.text(x0, 0.25, "\n".join(defs), fontsize=7, va="top", family="monospace")
        # a rebuild costs a fraction of a second and dragging a slider fires a
        # change on every mouse move, so changes only arm this timer and the
        # rebuild runs once, a moment after the slider stops moving
        self.timer = self.fig.canvas.new_timer(interval=300)
        self.timer.single_shot = True
        self.timer.add_callback(self.update, None)
        self.update(None)

    def schedule(self, _):
        self.status.set_text("...")
        self.timer.stop()
        self.timer.start()

    def values(self):
        v = {k: s.val for k, s in self.sliders.items()}
        for k in ("nrow", "ncol", "range", "seed", "n_points", "n_paths", "n_windows", "window", "num_reals"):
            v[k] = int(round(v[k]))
        return v

    def build(self):
        v = self.values()
        cfg = af.make_config(v["nrow"], v["ncol"], v["n_points"], v["n_paths"], v["n_windows"],
                             v["window"], v["a"], v["seed"])
        x = af.draw_fields(v["nrow"], v["ncol"], 1, v["range"], v["sill"], v["seed"])[0]
        return v, cfg, x

    def update(self, _):
        v, cfg, x = self.build()
        y = af.forward(x, cfg)
        npt, npa = len(cfg["points"]), len(cfg["paths"])

        ax = self.ax_field
        ax.clear()
        lim = max(1.0, np.abs(x).max())
        im = ax.imshow(x, cmap="viridis", vmin=-lim, vmax=lim, origin="upper", interpolation="nearest")
        if not hasattr(self, "cbar"):
            self.cbar = self.fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
            self.cbar.set_label("truth log10 k", fontsize=8)
        else:
            self.cbar.update_normal(im)
        for cells in cfg["paths"]:
            cells = np.asarray(cells)
            ax.plot(cells[:, 1], cells[:, 0], color="w", lw=1.2, alpha=0.9)
        w = cfg["window"]
        for i, j in cfg["windows"]:
            ax.add_patch(matplotlib.patches.Rectangle((j - 0.5, i - 0.5), w, w, fill=False,
                                                      ec="orange", lw=1.3))
        if npt:
            pts = np.asarray(cfg["points"])
            ax.scatter(pts[:, 1], pts[:, 0], s=22, c="r", edgecolor="k", lw=0.4, zorder=5)
        ax.set_title("truth field: red = point obs, white = paths, orange = blocks", fontsize=9)
        ax.set_xticks([])
        ax.set_yticks([])

        ax = self.ax_point
        ax.clear()
        xx = np.linspace(-3.0 * np.sqrt(v["sill"]), 3.0 * np.sqrt(v["sill"]), 200)
        ax.plot(xx, xx, color="0.6", lw=0.8, label="linear")
        ax.plot(xx, xx if v["a"] == 0 else np.sinh(v["a"] * xx) / v["a"], color="r", lw=1.2,
                label="sinh(ax)/a")
        if npt:
            pts = np.asarray(cfg["points"])
            ax.scatter(x[pts[:, 0], pts[:, 1]], y[:npt], s=14, c="r", edgecolor="k", lw=0.4, zorder=5)
        ax.set_xlabel("log10 k at the cell", fontsize=8)
        ax.set_ylabel("point obs", fontsize=8)
        ax.set_title("point obs vs the cell value", fontsize=9)
        ax.legend(fontsize=7, loc="upper left")
        ax.tick_params(labelsize=7)
        ax.grid(alpha=0.3)

        ax = self.ax_path
        ax.clear()
        lin, nl, kind = [], [], []
        for cells in cfg["paths"]:
            cells = np.asarray(cells)
            lin.append(x[cells[:, 0], cells[:, 1]].mean())
        nl += list(y[npt:npt + npa])
        kind += ["path"] * npa
        for i, j in cfg["windows"]:
            lin.append(x[i:i + w, j:j + w].mean())
        nl += list(y[npt + npa:])
        kind += ["block"] * len(cfg["windows"])
        lin, nl, kind = np.array(lin), np.array(nl), np.array(kind)
        if len(lin):
            r = max(np.abs(lin).max(), np.abs(nl).max()) * 1.1
            ax.plot([-r, r], [-r, r], color="0.6", lw=0.8)
            m = kind == "path"
            ax.scatter(lin[m], nl[m], s=16, c="w", edgecolor="k", lw=0.6, label="path (soft min, below the line)")
            ax.scatter(lin[~m], nl[~m], s=16, c="orange", edgecolor="k", lw=0.4, label="block (soft max, above)")
            ax.legend(fontsize=7, loc="upper left")
        ax.set_xlabel("mean log10 k of the cells", fontsize=8)
        ax.set_ylabel("obs", fontsize=8)
        ax.set_title("path and block obs vs the linear average", fontsize=9)
        ax.tick_params(labelsize=7)
        ax.grid(alpha=0.3)

        npar = v["nrow"] * v["ncol"]
        nobs = npt + npa + len(cfg["windows"])
        cov = "yes" if npar <= 5000 else "no (>5000 pars)"
        self.status.set_text("")
        self.summary.set_text("parameters  {0:>8d}\nobservations{1:>8d}\nprior cov   {2}\n"
                              "mean path   {3:5.1f} cells".format(
                                  npar, nobs, cov,
                                  np.mean([len(c) for c in cfg["paths"]]) if npa else 0.0))
        self.fig.canvas.draw_idle()

    def write(self, _):
        v = self.values()
        self.status.set_text("writing...")
        self.fig.canvas.draw()
        af.setup(new_d=self.out_d, nrow=v["nrow"], ncol=v["ncol"], n_points=v["n_points"],
                 n_paths=v["n_paths"], n_windows=v["n_windows"], window=v["window"], a=v["a"],
                 range_cells=float(v["range"]), sill=v["sill"], noise_std=v["noise"],
                 num_reals=v["num_reals"], seed=v["seed"])
        self.status.set_text("wrote ./{0}".format(self.out_d))
        self.fig.canvas.draw_idle()


if __name__ == "__main__":
    Dashboard()
    plt.show()
