import sys
sys.path.insert(0, '/Volumes/T7Shield/ncfa/athena/vis/python')
import athena_read

import matplotlib.pyplot as plt
import numpy as np

import data_handler


dir = "/Volumes/T7Shield/ncfa/athena/out/"
paths = []
for i in range(3):
    paths.append(dir + f"frb.out1.{i:05d}.athdf")
    data = athena_read.athdf(paths[i])

    x1v = np.log10(data['x1v'])
    x2v = np.log10(data['x2v'])
    x3v = np.log10(data['x3v'])
    Bcc1 = np.log10(data['Bcc1'])
    Bcc2 = np.log10(data['Bcc2'])
    Bcc3  = data['Bcc3'][0]
    _Bcc3 = np.copy(Bcc3) * -1
    Bcc3[Bcc3 < 0] = 1e-16
    Bcc3 = np.log10(Bcc3)
    
    
    _Bcc3[_Bcc3 < 0] = 1e-16
    _Bcc3 = np.log10(_Bcc3)

    fig, (ax1, ax2) = plt.subplots(1,2)
    im1 = ax1.imshow(Bcc3, interpolation='nearest', origin='lower', cmap='plasma', vmax=-5, vmin=-10)
    fig.colorbar(im1, ax=ax1, fraction=0.046, cmap='plasma', label=r'log$_{10}$B')
    ax1.set_xlabel('X'); ax1.set_ylabel('Y')
 
    im2 = ax2.imshow(_Bcc3, interpolation='nearest', origin='lower', cmap='plasma', vmax=-5, vmin=-10)
    fig.colorbar(im2, ax=ax2, fraction=0.046, cmap='plasma', label=r'-log$_{10}$B')
    ax2.set_xlabel('X'); ax2.set_ylabel('Y')

    plt.tight_layout()
    plt.savefig(f"out/{i:05d}.pdf", bbox_inches='tight')
    plt.close()