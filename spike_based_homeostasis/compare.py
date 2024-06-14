#!/usr/bin/env python3
"""Plots to compare results from Brian 2 and Arbor
"""

import json
import numpy as np
import matplotlib.pyplot as plt
import sklearn.metrics

def compute_and_print_goodness(data_1, data_2):
    """
    Computes and prints R^2 and RMSE measures for two given datasets (uses two different ways to compute each measure).
    """
    # compute R^2
    rsq = sklearn.metrics.r2_score(data_1, data_2)
    print(f"  R^2 = {rsq}")
    
    # compute RMSE
    rmse = np.sqrt(sklearn.metrics.mean_squared_error(data_1, data_2))
    print(f"  RMSE = {rmse}")
    

def plot(config,
         arbor_input, arbor_traces, arbor_spikes,
         brian2_input, brian2_traces, brian2_spikes):
    """
    Plots comparisons of spikes and traces between Brian2 and Arbor (using first two default colors)
    """

    fig, axes = plt.subplots(3, sharex=True)

    # plot comparison of input rate traces
    axes[0].plot(arbor_input[:, 0], arbor_input[:, 1], label='Arbor', color='#1f77b4')
    axes[0].plot(brian2_input[:, 0], brian2_input[:, 1], label='Brian', linestyle='dashed', color='#ff7f0e')
    axes[0].set_ylabel("Varying input rate (Hz)")

    # plot comparison of synaptic weight traces
    axes[1].plot(arbor_traces[:, 0], arbor_traces[:, 4] + config["stimulus"]["steady"]["weight"],
                 label='Arbor', color='#1f77b4')
    axes[1].plot(brian2_traces[:, 0], brian2_traces[:, 2] + config["stimulus"]["steady"]["weight"],
                 label='Brian', color='#ff7f0e', linestyle='dashed')
    #axes[1].set_ylabel("Synaptic weight (µS)")
    axes[1].set_ylabel("Synaptic weight (arb. units)")

    # compute actual rates and plot comparison
    dt = 1000 #100 # in ms
    times = np.arange(0, config["simulation"]["runtime"], dt)
    idxs = np.searchsorted(arbor_spikes, times)
    counts = np.array([len(arbor_spikes[start:stop]) for start, stop in zip(idxs, idxs[1:])])
    arbor_firing_rate = counts / (dt/1000) # convert to Hz
    idxs = np.searchsorted(brian2_spikes, times)
    counts = np.array([len(brian2_spikes[start:stop]) for start, stop in zip(idxs, idxs[1:])])
    brian2_firing_rate = counts / (dt/1000) # convert to Hz
    axes[2].plot(times[:-1], arbor_firing_rate, label='Arbor', color='#1f77b4')
    axes[2].plot(times[:-1], brian2_firing_rate, label='Brian', linestyle='dashed', color='#ff7f0e')
    axes[2].set_ylabel("Target neuron rate (Hz)")

    axes[2].set_xlabel("Time (ms)")

    fig.legend(*axes[2].get_legend_handles_labels(), loc="lower right")

    # compute and print R^2 and RMSE
    compute_and_print_goodness(arbor_traces[:, 4], brian2_traces[:, 2])

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

    arbor_input = np.loadtxt("arbor_input.dat")
    arbor_traces = np.loadtxt("arbor_traces.dat")
    arbor_spikes = np.loadtxt("arbor_spikes.dat")
    brian2_input = np.loadtxt("brian2_input.dat")
    brian2_traces = np.loadtxt("brian2_traces.dat")
    brian2_spikes = np.loadtxt("brian2_spikes.dat")
    config = json.load(open(args.config, 'r'))

    fig = plot(config, arbor_input, arbor_traces, arbor_spikes, brian2_input, brian2_traces, brian2_spikes)

    if args.save:
        fig.savefig(args.save)

    if args.show:
        plt.show()
