#!/usr/bin/env python3

import json
import numpy as np
import matplotlib.pyplot as plt

def main(variant):

    config = json.load(open(f"config_{variant}.json"))

    # um^2 -> cm^2
    surface = 4 * np.pi * config["neuron"]["radius"]**2 * 1e-6**2/1e-2**2

    data = np.loadtxt(f'arbor_traces_{variant}.dat')
    t = data[:, 0]
    vm = data[:, 1]

    # mA/cm^2 -> nA
    w = data[:, 2]*1e-3 * surface / 1e-9

    fig = plt.figure(figsize=(10, 5))
    fig.suptitle(config["pattern"])
    gs = fig.add_gridspec(2, 2)

    ax_vm = fig.add_subplot(gs[0, 0])
    ax_w = fig.add_subplot(gs[1, 0])
    ax_vm_w = fig.add_subplot(gs[:, 1])

    ax_vm.plot(t, vm)
    ax_w.plot(t, w)
    ax_vm_w.plot(vm, w)

    ax_vm.tick_params(labelbottom=False)

    ax_vm.set_ylabel("V [mV]")

    ax_w.set_xlabel("t [ms]")
    ax_w.set_ylabel("w [nA]")

    ax_vm_w.set_xlabel("V [mV]")
    ax_vm_w.set_ylabel("w [nA]")

    ax_vm_w.yaxis.tick_right()
    ax_vm_w.yaxis.set_label_position("right")

    fig.savefig(f"{variant}.png")

if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('variant', help="name of variant, e.g., tonic_spiking")
    args = parser.parse_args()

    main(args.variant)
