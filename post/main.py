import sys
sys.path.insert(0, '/Volumes/T7Shield/ncfa/athena/vis/python')
import athena_read

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy import interpolate
from scipy import integrate

import data_handler

# How to track NS position to draw los?

def main():
    dir = "./data/"
    paths = []

    for t in range(1):
        paths.append(dir + f"frb.out1.{t:05d}.athdf")
        raw_data = athena_read.athdf(paths[t])

        params = ['x1v', 'x2v', 'x3v', 'Bcc1', 'Bcc2', 'Bcc3', 'rho']
        data = parse_data(raw_data, params)
        # plot_b_vec(data)
        # plot_b_xyz(data, i)
        theta_los = np.pi/4
        phi_los = 0
        calc_B_parr(data, theta_los, phi_los, t)

def parse_data(raw_data, params):
    data = {}
    for param in params:
        data[param] = raw_data[param]
    data['x3v'] = np.array(np.size(data['x1v'])*[np.array(np.size(data['x1v'])*[0])]) # artifically insert
    data['x1v'], data['x2v'] = np.meshgrid(data['x1v'], data['x2v'])
    # print(data)
    return data

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

    x, y, z = np.meshgrid(x1v[:50],x2v[:50],x3v[:50])

    # vecs = [[[x1v[i][j], x2v[i][j], x3v[i][j], Bcc1[i][j], Bcc2[i][j], Bcc3[i][j]] for j,_ in enumerate(x1v[i])] for i,_ in enumerate(x1v)]
    # print(vecs)


    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.quiver(x, y, z, Bcc1, Bcc2, Bcc3)
    ax.set_xlim([-2,2])
    ax.set_ylim([-2, 2])
    # ax.set_zlim([-1, 8])
    plt.show()

def plot_b_xyz(data, i):
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

    vlim=[-8,-5]

    fig, axs = plt.subplots(3,2)
    im1 = axs[0,0].imshow(Bcc1, interpolation='nearest', origin='lower', cmap='plasma', vmax=vlim[0], vmin=vlim[1])
    fig.colorbar(im1, ax=axs[0,0], fraction=0.046, label=r'$\log_{10}B_x$')
    axs[0,0].set_xlabel('X'); axs[0,0].set_ylabel('Y')

    im2 = axs[0,1].imshow(_Bcc1, interpolation='nearest', origin='lower', cmap='plasma', vmax=vlim[0], vmin=vlim[1])
    fig.colorbar(im2, ax=axs[0,1], fraction=0.046, label=r'-$\log_{10}B_x$')
    axs[0,1].set_xlabel('X'); axs[0,1].set_ylabel('Y')

    im3 = axs[1,0].imshow(Bcc2, interpolation='nearest', origin='lower', cmap='plasma', vmax=vlim[0], vmin=vlim[1])
    fig.colorbar(im3, ax=axs[1,0], fraction=0.046, label=r'$\log_{10}B_y$')
    axs[1,0].set_xlabel('X'); axs[1,0].set_ylabel('Y')

    im4 = axs[1,1].imshow(_Bcc2, interpolation='nearest', origin='lower', cmap='plasma', vmax=vlim[0], vmin=vlim[1])
    fig.colorbar(im4, ax=axs[1,1], fraction=0.046, label=r'-$\log_{10}B_y$')
    axs[1,1].set_xlabel('X'); axs[1,1].set_ylabel('Y')

    im5 = axs[2,0].imshow(Bcc3, interpolation='nearest', origin='lower', cmap='plasma', vmax=vlim[0], vmin=vlim[1])
    fig.colorbar(im5, ax=axs[2,0], fraction=0.046, label=r'$\log_{10}B_z$')
    axs[2,0].set_xlabel('X'); axs[2,0].set_ylabel('Y')

    im6 = axs[2,1].imshow(_Bcc3, interpolation='nearest', origin='lower', cmap='plasma', vmax=vlim[0], vmin=vlim[1])
    fig.colorbar(im6, ax=axs[2,1], fraction=0.046, label=r'-$\log_{10}B_z$')
    axs[2,1].set_xlabel('X'); axs[2,1].set_ylabel('Y')

    # plt.suptitle('Magnetic Field in X-Y-Z', x=(fig.subplotpars.right + fig.subplotpars.left)/2)

    plt.tight_layout()
    plt.savefig(f"out/{i:05d}.pdf", bbox_inches='tight')
    plt.close()

def calc_B_parr(data, theta_los, phi_los, t):
    # Setup numerical parameters
    l = 2
    n = 50
    dl = l/n
    
    # Initialize interpolation functions
    interp_B_x = interpolate.interp2d(data['x1v'][0], np.array([row[0] for row in data['x2v']]), data['Bcc1'][0], kind='cubic')
    interp_B_y = interpolate.interp2d(data['x1v'][0], np.array([row[0] for row in data['x2v']]), data['Bcc2'][0], kind='cubic')
    interp_B_z = interpolate.interp2d(data['x1v'][0], np.array([row[0] for row in data['x2v']]), data['Bcc3'][0], kind='cubic')
    interp_rho = interpolate.interp2d(data['x1v'][0], np.array([row[0] for row in data['x2v']]), data['rho'][0], kind='cubic')
    
    # Calculate grid spacing
    delta_x = dl * np.cos(theta_los)
    delta_y = dl * np.sin(theta_los)
    
    # Initialize arrays to store x and y coordinates
    x = np.zeros(n*2+1)
    y = np.zeros(n*2+1)

    # los unit vector
    n_los = [delta_x, delta_y] / np.sqrt(delta_x * delta_x + delta_y * delta_y)
    
    # Initialize los data arrays
    los_B_x = []
    los_B_y = []
    los_B_z = []
    los_rho = []

    # Iterate over los grid
    for i in range(-n, n+1):
        # Get los coordinates
        x[-(i+n)-1] = i * delta_x
        y[-(i+n)-1] = i * delta_y

        # Interpolate data along los
        los_B_x.append(interp_B_x(x[-(i+n)-1], y[-(i+n)-1])[0])
        los_B_y.append(interp_B_y(x[-(i+n)-1], y[-(i+n)-1])[0])
        los_B_z.append(interp_B_z(x[-(i+n)-1], y[-(i+n)-1])[0])
        los_rho.append(interp_rho(x[-(i+n)-1], y[-(i+n)-1])[0])

    # Calculate B_parr
    B_parr = 0
    for i in range(len(los_B_x)):
        B_x = los_B_x[i]
        B_y = los_B_y[i]
        rho = los_rho[i]
        # By dot product of los unit vector and B vector
        B_parr += np.dot(n_los, [B_x, B_y])

    # Normalize
    B_parr /= n
    print(B_parr)


    # plots

    fig, ax = plt.subplots(2,2, figsize=(10,8))
    ax[0,0].set_xlabel('X'); ax[0,0].set_ylabel('Y')
    ax[0,0].title.set_text('B_x')
    ax[0,0].set_xticks(np.arange(-2, 2.1, step=1)); ax[0,0].set_yticks(np.arange(-2, 2.1, step=1))
    im1 = ax[0,0].imshow(data['Bcc1'][0], origin='lower', extent=[-2,2,-2,2], cmap='coolwarm', vmax=1e-5, vmin=-1e-5)
    ax[0,0].arrow(x[1], y[1], (x[-1] - x[0])/10, (y[-1] - y[0])/10, head_width=0.1, color='grey')
    ax[0,0].plot(x, y, color='grey', alpha=0.5, ls='--')
    fig.colorbar(im1, ax=ax[0,0], fraction=0.046, label=r'$B_x$')

    ax[0,1].set_xlabel('X'); ax[0,1].set_ylabel('Y')
    ax[0,1].title.set_text('B_y')
    ax[0,1].set_xticks(np.arange(-2, 2.1, step=1)); ax[0,1].set_yticks(np.arange(-2, 2.1, step=1))
    im2 = ax[0,1].imshow(data['Bcc2'][0], origin='lower', extent=[-2,2,-2,2], cmap='coolwarm', vmax=1e-5, vmin=-1e-5)
    ax[0,1].arrow(x[1], y[1], (x[-1] - x[0])/10, (y[-1] - y[0])/10, head_width=0.1, color='grey')
    ax[0,1].plot(x, y, color='grey', alpha=0.5, ls='--')
    fig.colorbar(im2, ax=ax[0,1], fraction=0.046, label=r'$B_y$')

    ax[1,0].set_xlabel('X'); ax[1,0].set_ylabel('Y')
    ax[1,0].title.set_text('B_z')
    ax[1,0].set_xticks(np.arange(-2, 2.1, step=1)); ax[1,0].set_yticks(np.arange(-2, 2.1, step=1))
    im3 = ax[1,0].imshow(data['Bcc3'][0], origin='lower', extent=[-2,2,-2,2], cmap='coolwarm', vmax=1e-5, vmin=-1e-5)
    ax[1,0].arrow(x[1], y[1], (x[-1] - x[0])/10, (y[-1] - y[0])/10, head_width=0.1, color='grey')
    ax[1,0].plot(x, y, color='grey', alpha=0.5, ls='--')
    fig.colorbar(im3, ax=ax[1,0], fraction=0.046, label=r'$B_z$')

    ax[1,1].set_xlabel('X'); ax[1,1].set_ylabel('Y')
    ax[1,1].title.set_text(r'$\rho$')
    ax[1,1].set_xticks(np.arange(-2, 2.1, step=1)); ax[1,1].set_yticks(np.arange(-2, 2.1, step=1))
    im4 = ax[1,1].imshow(np.log10(data['rho'][0]), origin='lower', extent=[-2,2,-2,2], cmap='coolwarm')
    ax[1,1].arrow(x[1], y[1], (x[-1] - x[0])/10, (y[-1] - y[0])/10, head_width=0.1, color='grey')
    ax[1,1].plot(x, y, color='grey', alpha=0.5, ls='--')
    fig.colorbar(im4, ax=ax[1,1], fraction=0.046, label=r'$\log_{10}rho$')

    plt.tight_layout()
    plt.savefig(f"out/{t:05d}.pdf", bbox_inches='tight')
    plt.show()

    return B_parr



if __name__ == '__main__':
    main()