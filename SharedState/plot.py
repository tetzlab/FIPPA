#!/usr/bin/env python3

import json
import numpy as np
import matplotlib.pyplot as plt

def main():

    data = np.loadtxt(f'arbor_traces.dat')
    t = data[:, 0]
    vm = data[:, 1]
    pi = data[:, 2]

    fig, axes = plt.subplots(2)

    axes[0].plot(t, vm)
    axes[1].plot(t, pi)

    axes[0].set_ylabel("V [mV]")
    axes[1].set_ylabel("pi")

    axes[1].set_xlabel("t [ms]")

    fig.savefig(f"plot.png")

if __name__ == '__main__':

    main()
