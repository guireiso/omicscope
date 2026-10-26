"""OmicScope 2.0 logo: flat redraw of the original (nebula disc, orbits with nodes, telescope) + wordmark."""
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, Ellipse, FancyBboxPatch, Polygon
from matplotlib.transforms import Affine2D

plt.rcParams['svg.fonttype'] = 'path'
OUT = sys.argv[1]
NAVY, DEEP, TEAL = '#14304f', '#0b1a2e', '#1f6f8b'
PEACH, CORAL, INK = '#efc69b', '#e76f51', '#0b1a2e'


def disc(ax):
    # background: radial gradient navy -> teal ring -> warm core, clipped to the disc
    n = 600
    y, x = np.mgrid[-1:1:n * 1j, -1:1:n * 1j]
    r = np.hypot(x - 0.05, y + 0.02)
    cols = np.array([[0.91, 0.45, 0.36], [0.55, 0.30, 0.42], [0.12, 0.36, 0.50], [0.08, 0.19, 0.31], [0.05, 0.11, 0.19]])
    stops = np.array([0.0, 0.28, 0.55, 0.8, 1.0])
    img = np.stack([np.interp(r, stops, cols[:, c]) for c in range(3)], axis=-1)
    im = ax.imshow(img, extent=(-1, 1, -1, 1), origin='lower', zorder=1, interpolation='bicubic')
    clip = Circle((0, 0), 1, transform=ax.transData)
    im.set_clip_path(clip)
    ax.add_patch(Circle((0, 0), 1, fill=False, lw=6, ec=DEEP, zorder=6))
    # a few stars, kept sparse
    rng = np.random.default_rng(3)
    a, rr = rng.uniform(0, 2 * np.pi, 26), np.sqrt(rng.uniform(0.35, 0.92, 26))
    ax.scatter(rr * np.cos(a), rr * np.sin(a), s=rng.uniform(1.5, 7, 26), c='#cfe8f3', alpha=0.85, lw=0, zorder=2)
    # orbits
    for ang, w, h in [(25, 1.38, 0.62), (-35, 1.32, 0.58), (90, 1.2, 0.52)]:
        ax.add_patch(Ellipse((0.04, 0.0), w, h, angle=ang, fill=False, ec=PEACH, lw=2.2, alpha=0.9, zorder=3))
    # network nodes on the orbits, a few edges
    def on(ang, w, h, t):
        p = np.array([w / 2 * np.cos(t), h / 2 * np.sin(t)])
        c, s = np.cos(np.radians(ang)), np.sin(np.radians(ang))
        return np.array([c * p[0] - s * p[1] + 0.04, s * p[0] + c * p[1]])
    nodes = [on(25, 1.38, 0.62, 0.35), on(-35, 1.32, 0.58, 2.6), on(90, 1.2, 0.52, 1.2),
             on(25, 1.38, 0.62, 3.6), on(-35, 1.32, 0.58, 5.5)]
    for i, j in [(0, 2), (2, 1), (0, 4)]:
        ax.plot(*zip(nodes[i], nodes[j]), color=PEACH, lw=1.4, alpha=0.75, zorder=3)
    for k, (px, py) in enumerate(nodes):
        ax.add_patch(Circle((px, py), 0.052 if k != 3 else 0.042, fc=CORAL if k == 3 else PEACH, ec='none', zorder=4))
    # telescope: tube, eyepiece and tripod, as one dark silhouette
    tube = Affine2D().rotate_deg(38).translate(0.02, 0.18)
    t = Polygon([(-0.5, -0.05), (0.42, -0.075), (0.42, 0.075), (-0.5, 0.05)], closed=True)
    t.set_transform(tube + ax.transData)
    ax.add_patch(t)
    t.set(fc=INK, ec='none', zorder=5)
    cap = Polygon([(0.42, -0.092), (0.49, -0.092), (0.49, 0.092), (0.42, 0.092)], closed=True, fc=INK, ec='none', zorder=5)
    cap.set_transform(tube + ax.transData)
    ax.add_patch(cap)
    eye = Polygon([(-0.6, -0.032), (-0.5, -0.036), (-0.5, 0.036), (-0.6, 0.032)], closed=True, fc=INK, ec='none', zorder=5)
    eye.set_transform(tube + ax.transData)
    ax.add_patch(eye)
    hx, hy = 0.0, 0.1
    ax.add_patch(Circle((hx, hy), 0.045, fc=INK, ec='none', zorder=5))
    for fx in (-0.55, 0.0, 0.55):
        ax.plot([hx, fx], [hy, -0.9 if fx == 0 else -0.86], color=INK, lw=4.5, solid_capstyle='round', zorder=5)
    ax.plot([-0.33, 0.33], [-0.38, -0.38], color=INK, lw=2.4, zorder=5)


def _bbox(ax, artist):
    fig = ax.figure
    fig.canvas.draw()
    return artist.get_window_extent(fig.canvas.get_renderer()).transformed(ax.transData.inverted())


def wordmark(ax, x, y, size, ha='center', color=INK):
    from matplotlib import patheffects
    return ax.text(x, y, 'OmicScope', fontsize=size, fontfamily='Bahnschrift', color=color, ha=ha, va='center',
                   zorder=7, path_effects=[patheffects.withStroke(linewidth=size / 22, foreground=color)])


def badge(ax, word, size, gap=0.05, rise=0.45):
    """'2.0' pill just after the wordmark, raised like an exponent; sized from its own text."""
    wb = _bbox(ax, word)
    t = ax.text(0, 0, '2.0', fontsize=size, fontfamily='Bahnschrift', color='white', ha='left', va='center', zorder=8)
    tb = _bbox(ax, t)
    w, h = tb.width, tb.height
    padx, pady = 0.45 * h, 0.32 * h
    x0 = wb.x1 + gap + padx
    yc = wb.y0 + wb.height * (0.5 + rise)
    t.set_position((x0, yc))
    ax.add_patch(FancyBboxPatch((x0 - padx, yc - h / 2 - pady), w + 2 * padx, h + 2 * pady,
                                boxstyle=f'round,pad=0,rounding_size={(h + 2 * pady) / 2}', fc=CORAL, ec='none',
                                zorder=7))


def figure(w, h, xlim, ylim):
    fig = plt.figure(figsize=(w, h))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect('equal')
    ax.axis('off')
    return fig, ax


# stacked logo (disc + wordmark + badge)
fig, ax = figure(6, 7.0, (-1.15, 1.15), (-1.62, 1.08))
disc(ax)
badge(ax, wordmark(ax, -0.12, -1.3, 50), 20, gap=0.04, rise=0.32)
for ext in ('svg', 'png'):
    fig.savefig(f'{OUT}/omicscope_2_logo.{ext}', dpi=300, transparent=True)
plt.close(fig)

# icon (disc + badge), for favicons and the sidebar
fig, ax = figure(4, 4, (-1.08, 1.08), (-1.08, 1.08))
disc(ax)
t = ax.text(0.62, -0.86, '2.0', fontsize=20, fontfamily='Bahnschrift', color='white', ha='center', va='center',
            zorder=8)
tb = _bbox(ax, t)
ax.add_patch(FancyBboxPatch((tb.x0 - 0.45 * tb.height, tb.y0 - 0.32 * tb.height), tb.width + 0.9 * tb.height,
                            tb.height * 1.64, boxstyle=f'round,pad=0,rounding_size={tb.height * 0.82}', fc=CORAL,
                            ec='none', zorder=7))
for ext in ('svg', 'png'):
    fig.savefig(f'{OUT}/omicscope_2_icon.{ext}', dpi=200, transparent=True)
plt.close(fig)

# horizontal lock-up for headers
fig, ax = figure(8.6, 2.6, (-1.1, 6.1), (-1.08, 1.08))
disc(ax)
badge(ax, wordmark(ax, 1.3, -0.04, 60, ha='left'), 22, gap=0.08, rise=0.3)
for ext in ('svg', 'png'):
    fig.savefig(f'{OUT}/omicscope_2_horizontal.{ext}', dpi=200, transparent=True)
plt.close(fig)
print('ok')
