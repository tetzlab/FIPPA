#!/usr/bin/env python3

import os
import json
import numpy as np
import matplotlib.pyplot as plt

def main():

    data = np.loadtxt(f'arbor_traces_bac.dat')
    t = data[:, 0]
    vm_soma = data[:, 1]
    vm_proximal = data[:, 2]
    vm_distal = data[:, 3]

    """
    fig, ax = plt.subplots(3, figsize=(10, 5))

    ax[0].plot(t, vm_soma)
    ax[1].plot(t, vm_proximal)
    ax[2].plot(t, vm_distal)
    """

    fig, ax = plt.subplots(1, figsize=(10, 5))

    ax.plot(t, vm_soma, label="soma", color='tab:blue')
    ax.plot(t, vm_proximal, label="proximal", color='tab:orange')
    ax.plot(t, vm_distal, label="distal", color='tab:gray')

    plt.legend()

    ax.set_xlabel("t [ms]")
    ax.set_ylabel("V [mV]")

    #ax.set_ylim(-100, 100)

    dir_path = os.getcwd().split('/')[-1]

    fig.savefig(f"{dir_path}.png")

if __name__ == '__main__':

    main()
