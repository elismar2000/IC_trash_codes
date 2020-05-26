import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

plt.rc('grid', color='#397939', linewidth=1, linestyle='-')
plt.rc('xtick', labelsize=10)
plt.rc('ytick', labelsize=5)

width, height = matplotlib.rcParams['figure.figsize']
size = min(width, height)
fig = plt.figure(figsize=(size, size))
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], polar=True, facecolor='#cfd98c')

ax.set_rmax(20.0)
plt.grid(True)

title = ax.text(0.5,0.85, "", bbox={'facecolor':'w', 'alpha':0.5, 'pad':5},
                transform=ax.transAxes, ha="center")

def data_gen(t=0):
    tw = 10
    phase = 0
    while True:
        if phase < 2*180:
            phase += 2
        else:
            phase=0
        yield tw, phase

def update(data):
    tw, phase = data
    title.set_text(u"|TW| = {}, Angle: {}°".format(tw, phase))
    arr1 = ax.arrow(np.deg2rad(phase), 0, 0, tw, alpha = 0.5, width = 0.080,
             edgecolor = 'red', facecolor = 'red', lw = 2, zorder = 5)
    return arr1,title,

ani = animation.FuncAnimation(fig, update, data_gen, interval=100, blit=True, repeat=False)

plt.show()
