#!/usr/bin/env python3

import json

import numpy as np
import matplotlib.pyplot as plt


def plot(config, traces, spikes):

    fig, axes = plt.subplots(2, sharex=True)

    dt = 100 # in ms
    times = np.arange(0, config["simulation"]["runtime"], dt)
    idxs = np.searchsorted(spikes, times)
    counts = np.array([len(spikes[start:stop]) for start, stop in zip(idxs, idxs[1:])])
    firing_rate = counts / (dt/1000) # convert to Hz

    axes[0].plot(times[:-1], firing_rate)
    axes[0].set_ylabel("rate [Hz]")
    axes[1].plot(traces[:,0], traces[:,4] + config["stimulus"]["steady"]["weight"])
    axes[1].set_ylabel("weight [$\mu$S]")

    axes[1].set_xlabel("time [ms]")

    return fig

if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('config', help="name of config file")
    parser.add_argument('--show', help="show plot",
                        action="store_true", default=False)
    parser.add_argument('--save', help="save to given file name", default=None)

    # parse the command line arguments
    args = parser.parse_args()

    if not args.show and not args.save:
        print("Neither --show nor --save selected, "
              "simulation will run but no output will be produced.")

    spikes = np.loadtxt("spikes.dat")
    traces = np.loadtxt("traces.dat")
    config = json.load(open(args.config, 'r'))

    fig = plot(config, traces, spikes)

    if args.save:
        fig.savefig(args.save)

    if args.show:
        plt.show()
