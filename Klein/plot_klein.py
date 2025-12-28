import numpy as np
import glob, re, os
import matplotlib.pyplot as plt
import imageio

files = sorted(glob.glob("data/klein_*.csv"),
               key=lambda s: int(re.findall(r"(\d+)", s)[-1]))

if not files:
    raise RuntimeError("No data/klein_*.csv found. Run the C++ program first.")

data0 = np.loadtxt(files[0], delimiter=",")
x0 = data0[:,0]
total0 = data0[:,7]  # last column = total
xmin, xmax = x0.min(), x0.max()

step_x = xmin + 0.55*(xmax-xmin)

out_gif = "klein_paradox.gif"
fps = 25

tmp_pngs = []
for idx, f in enumerate(files):
    d = np.loadtxt(f, delimiter=",")
    x = d[:,0]
    p1 = d[:,5]
    p2 = d[:,6]
    total = d[:,7]

    fig, ax = plt.subplots(figsize=(9,4))
    ax.plot(x, total, lw=2, label=r"$|\psi|^2$")
    ax.plot(x, p1, alpha=0.7, label=r"$|\psi_1|^2$")
    ax.plot(x, p2, alpha=0.7, label=r"$|\psi_2|^2$")

    ax.axvspan(step_x, xmax, alpha=0.15, label="Step region (V0)")

    ax.axvline(step_x, lw=2)

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(0, max(total0.max(), total.max()) * 1.25)
    ax.set_xlabel("x")
    ax.set_ylabel("Probability density")
    ax.set_title("1D Dirac wavepacket at supercritical step (Klein tunneling)")
    ax.legend(loc="upper right")

    png = f"_tmp_{idx:04d}.png"
    fig.tight_layout()
    fig.savefig(png, dpi=140)
    plt.close(fig)
    tmp_pngs.append(png)

print("Building GIF...")
images = [imageio.imread(p) for p in tmp_pngs]
imageio.mimsave(out_gif, images, fps=fps)

for p in tmp_pngs:
    os.remove(p)

print("Saved:", out_gif)
