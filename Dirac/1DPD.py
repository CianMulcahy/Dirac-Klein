import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import re
from matplotlib.animation import FuncAnimation

data_dir = "data"        
pattern = os.path.join(data_dir, "psi_*.csv")
save_gif = True             
gif_filename = "dirac_wave.gif"
fps = 30            
skip = 1                

files = sorted(glob.glob(pattern), key=lambda f: int(re.findall(r"(\d+)", f)[-1]))

if not files:
    raise FileNotFoundError("No CSV files found in ./data/ — did you run the simulation first?")

x, psi1, psi2, total = np.loadtxt(files[0], delimiter=",", unpack=True)
N = len(x)

fig, ax = plt.subplots(figsize=(8, 5))
line_total, = ax.plot([], [], color="black", lw=2, label=r"$|\psi|^2$")
line_p1, = ax.plot([], [], color="blue", alpha=0.7, label=r"$|\psi_1|^2$")
line_p2, = ax.plot([], [], color="red", alpha=0.7, label=r"$|\psi_2|^2$")
ax.set_xlim(x.min(), x.max())
ax.set_ylim(0, max(total) * 1.5)
ax.set_xlabel("x")
ax.set_ylabel("Probability density")
ax.legend(loc="upper right")
ax.set_title("1D Dirac Wavepacket Evolution")

def update(frame):
    file = files[frame]
    x, psi1, psi2, total = np.loadtxt(file, delimiter=",", unpack=True)
    line_total.set_data(x, total)
    line_p1.set_data(x, psi1)
    line_p2.set_data(x, psi2)
    ax.set_ylim(0, max(total) * 1.5)
    ax.set_title(f"1D Dirac Wavepacket Evolution (t = {frame})")
    return line_total, line_p1, line_p2

ani = FuncAnimation(fig, update, frames=range(0, len(files), skip),
                    interval=1000/fps, blit=True)

if save_gif:
    print("Saving animation...")
    ani.save(gif_filename, fps=fps, writer="pillow")
    print(f"Saved animation as {gif_filename}")
else:
    plt.show()
