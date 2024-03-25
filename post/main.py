import matplotlib.pyplot as plt
import numpy as np

import data_handler


dir = "/Volumes/T7Shield/ncfa/athena/out/"
paths = []
for i in range(350):
    paths.append(dir + f"frb.out1.{i:05d}.athdf")
    data = data_handler.read(paths[i])
    data = np.log10(data)

    plt.figure()
    plt.imshow(data, interpolation='nearest', origin='lower', cmap='plasma', vmin=-10, vmax=-5)
    plt.colorbar(cmap='plasma', label=r'log$_{10}$B')
    plt.xlabel('X'); plt.ylabel('Y')
 
    plt.savefig(f"out/{i:05d}.png")
    plt.close()