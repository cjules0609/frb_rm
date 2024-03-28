import sys
sys.path.insert(0, '/Volumes/T7Shield/ncfa/athena/vis/python')
import athena_read

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

import data_handler


def main():
    dir = "/Volumes/T7Shield/ncfa/athena/out/"
    paths = []

    for i in range(2):
        paths.append(dir + f"frb.out1.{i:05d}.athdf")
        data = athena_read.athdf(paths[i])

        plot_b_vec(data)


def plot_b_vec(data):
    x1v = data['x1v']
    x1v = np.array(512*[x1v])
    x2v = data['x2v']
    x2v = np.array(512*[x2v])
    x3v = data['x3v']
    x3v = np.array(512*[np.array(512*[0])])
    Bcc1  = data['Bcc1'][0]
    Bcc2  = data['Bcc2'][0]
    Bcc3  = data['Bcc3'][0]

    # vecs = [[[x1v[i][j], x2v[i][j], x3v[i][j], Bcc1[i][j], Bcc2[i][j], Bcc3[i][j]] for j,_ in enumerate(x1v[i])] for i,_ in enumerate(x1v)]
    # print(vecs)


    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.quiver(x1v, x2v, x3v, Bcc1, Bcc2, Bcc3)
    # ax.set_xlim([-1, 0.5])
    # ax.set_ylim([-1, 1.5])
    # ax.set_zlim([-1, 8])
    plt.show()

def plot_b_xyz(data):
    Bcc1  = data['Bcc1'][0] + 1e-32
    _Bcc1 = np.copy(Bcc1) * -1
    Bcc1[Bcc1 < 0] = 1e-16
    Bcc1 = np.log10(Bcc1)
    _Bcc1[_Bcc1 < 0] = 1e-16
    _Bcc1 = np.log10(_Bcc1)

    Bcc2  = data['Bcc2'][0] + 1e-32
    _Bcc2 = np.copy(Bcc2) * -1
    Bcc2[Bcc2 < 0] = 1e-16
    Bcc2 = np.log10(Bcc2)
    _Bcc2[_Bcc2 < 0] = 1e-16
    _Bcc2 = np.log10(_Bcc2)

    Bcc3  = data['Bcc3'][0]
    _Bcc3 = np.copy(Bcc3) * -1
    Bcc3[Bcc3 < 0] = 1e-16
    Bcc3 = np.log10(Bcc3)
    _Bcc3[_Bcc3 < 0] = 1e-16
    _Bcc3 = np.log10(_Bcc3)

    fig, axs = plt.subplots(3,2)
    im1 = axs[0,0].imshow(Bcc1, interpolation='nearest', origin='lower', cmap='plasma', vmax=-5, vmin=-15)
    fig.colorbar(im1, ax=axs[0,0], fraction=0.046, label=r'$\log_{10}B_x$')
    axs[0,0].set_xlabel('X'); axs[0,0].set_ylabel('Y')

    im2 = axs[0,1].imshow(_Bcc1, interpolation='nearest', origin='lower', cmap='plasma', vmax=-5, vmin=-15)
    fig.colorbar(im2, ax=axs[0,1], fraction=0.046, label=r'-$\log_{10}B_x$')
    axs[0,1].set_xlabel('X'); axs[0,1].set_ylabel('Y')

    im3 = axs[1,0].imshow(Bcc2, interpolation='nearest', origin='lower', cmap='plasma', vmax=-5, vmin=-15)
    fig.colorbar(im3, ax=axs[1,0], fraction=0.046, label=r'$\log_{10}B_y$')
    axs[1,0].set_xlabel('X'); axs[1,0].set_ylabel('Y')

    im4 = axs[1,1].imshow(_Bcc2, interpolation='nearest', origin='lower', cmap='plasma', vmax=-5, vmin=-15)
    fig.colorbar(im4, ax=axs[1,1], fraction=0.046, label=r'-$\log_{10}B_y$')
    axs[1,1].set_xlabel('X'); axs[1,1].set_ylabel('Y')

    im5 = axs[2,0].imshow(Bcc3, interpolation='nearest', origin='lower', cmap='plasma', vmax=-5, vmin=-15)
    fig.colorbar(im5, ax=axs[2,0], fraction=0.046, label=r'$\log_{10}B_z$')
    axs[2,0].set_xlabel('X'); axs[2,0].set_ylabel('Y')

    im6 = axs[2,1].imshow(_Bcc3, interpolation='nearest', origin='lower', cmap='plasma', vmax=-5, vmin=-15)
    fig.colorbar(im6, ax=axs[2,1], fraction=0.046, label=r'-$\log_{10}B_z$')
    axs[2,1].set_xlabel('X'); axs[2,1].set_ylabel('Y')

    # plt.suptitle('Magnetic Field in X-Y-Z', x=(fig.subplotpars.right + fig.subplotpars.left)/2)

    plt.tight_layout()
    plt.savefig(f"out/{i:05d}.pdf", bbox_inches='tight')
    plt.close()

if __name__ == '__main__':
    main()