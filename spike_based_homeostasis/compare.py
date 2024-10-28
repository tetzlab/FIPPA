#!/usr/bin/env python3
"""Plots to compare results from Brian 2 and Arbor
"""

import json
import numpy as np
import matplotlib.pyplot as plt
import sklearn.metrics

def compute_and_print_goodness(data_1, data_2):
    """
    Computes and prints R^2, CV, and RMSE measures for two given datasets.
    """
    # compute R^2
    rsq = sklearn.metrics.r2_score(data_1, data_2)
    print(f"  R^2 = {rsq}")

	# compute coefficient of variation
    cv = sklearn.metrics.explained_variance_score(data_1, data_2)
    print(f"  CV = {cv}")
    
    # compute RMSE
    rmse = np.sqrt(sklearn.metrics.root_mean_squared_error(data_1, data_2))
    print(f"  RMSE = {rmse}")


def plot(config, num_trials,
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
    axes[1].plot(arbor_traces[:, 0], arbor_traces[:, 4],
                 label='Arbor', color='#1f77b4')
    axes[1].plot(brian2_traces[:, 0], brian2_traces[:, 2],
                 label='Brian', color='#ff7f0e', linestyle='dashed')
    #axes[1].set_ylabel("Synaptic weight (µS)")
    axes[1].set_ylabel("Synaptic weight (arb. units)")

    # compute actual rates and plot comparison
    time_window = 100 # ms
    time_bins = np.arange(0, config["simulation"]["runtime"], time_window)
    arbor_counts = np.array([])
    brian2_counts = np.array([])
    for time_bin in time_bins:
        current_count = len(arbor_spikes[np.logical_and(arbor_spikes >= time_bin, arbor_spikes < time_bin + time_window)])
        arbor_counts = np.append(arbor_counts, current_count)
        current_count = len(brian2_spikes[np.logical_and(brian2_spikes >= time_bin, brian2_spikes < time_bin + time_window)])
        brian2_counts = np.append(brian2_counts, current_count)
    #idxs = np.searchsorted(arbor_spikes, time_bins)
    #counts = np.array([len(arbor_spikes[start:stop]) for start, stop in zip(idxs, idxs[1:])])
    arbor_firing_rate = arbor_counts / num_trials / (time_window/1000) # convert to Hz
    #idxs = np.searchsorted(brian2_spikes, time_bins)
    #counts = np.array([len(brian2_spikes[start:stop]) for start, stop in zip(idxs, idxs[1:])])
    brian2_firing_rate = brian2_counts / num_trials / (time_window/1000) # convert to Hz
    axes[2].plot(time_bins, arbor_firing_rate, label='Arbor', color='#1f77b4')
    axes[2].plot(time_bins, brian2_firing_rate, label='Brian', linestyle='dashed', color='#ff7f0e')
    axes[2].set_ylim(-5, 105)
    axes[2].set_yticks(np.arange(0, 105, 25))
    axes[2].set_ylabel("Resulting rate (Hz)")

    axes[2].set_xlabel("Time (ms)")

    fig.legend(*axes[2].get_legend_handles_labels(), loc="lower right")

    # compute and print R^2 and RMSE for weight
    print("Goodness of fit (weight):")
    compute_and_print_goodness(arbor_traces[:, 4], brian2_traces[:, 2])

    # compute and print R^2 and RMSE for resulting firing rate
    print("Goodness of fit (resulting rate):")
    compute_and_print_goodness(arbor_firing_rate, brian2_firing_rate)

    return fig


if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('config', help="name of config file")
    parser.add_argument('num_trials', type=int, help="number of trials in the data")
    parser.add_argument('considered_neuron', type=int, help="the type of neuron to consider")
    parser.add_argument('--show', help="show plot",
                        action="store_true", default=False)
    parser.add_argument('--save', help="save to given file name", default=None)

    # parse the command line arguments
    args = parser.parse_args()

    if not args.show and not args.save:
        print("Neither --show nor --save selected, "
              "simulation will run but no output will be produced.")

    arbor_input = np.loadtxt(f"arbor_input_{args.considered_neuron}.dat")
    arbor_traces = np.loadtxt(f"arbor_traces_{args.considered_neuron}.dat")
    arbor_spikes = np.loadtxt(f"arbor_spikes_{args.considered_neuron}.dat")
    brian2_input = np.loadtxt(f"brian2_input_{args.considered_neuron}.dat")
    brian2_traces = np.loadtxt(f"brian2_traces_{args.considered_neuron}.dat")
    brian2_spikes = np.loadtxt(f"brian2_spikes_{args.considered_neuron}.dat")
    config = json.load(open(args.config, 'r'))

    fig = plot(config, args.num_trials, arbor_input, arbor_traces, arbor_spikes, brian2_input, brian2_traces, brian2_spikes)

    if args.save:
        fig.savefig(args.save)

    if args.show:
        plt.show()
