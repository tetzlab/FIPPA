#!/usr/bin/env python3
"""Plots to compare results from Brian 2 and Arbor
"""

import json
import numpy as np
import matplotlib.pyplot as plt


def main(variant):
    """Plots comparisons of spikes and traces between Brian2 and Arbor
    """
    config_lif = json.load(open(f"config_{variant}_lif.json"))
    config_classical = json.load(open(f"config_{variant}_classical.json"))

    fig, axes = plt.subplots(nrows=3, ncols=2, sharex=False, figsize=(10, 10))

    stimulus_times_exc = config_lif["synapses"]["cond_exp_stdp"]["stimulus_times"]
    stimulus_times_inh = config_lif["synapses"]["cond_exp"]["stimulus_times"]

    # plot trace comparisons for LIF neurons
    arbor_data = np.loadtxt(f'arbor_traces_{variant}_lif.dat')
    brian_data = np.loadtxt(f'brian2_traces_{variant}_lif.dat')

    # time, membrane voltage, excitatory conductance, inhibitory conductance,
    # weight
    for i in range(1, 5):

        for spike_time in stimulus_times_exc:
            axes.flat[i - 1].axvline(spike_time, ymin=0.9, ymax=1, color='blue')

        for spike_time in stimulus_times_inh:
            axes.flat[i - 1].axvline(spike_time, ymin=0.85,
                                     ymax=0.95, color='red')

        # offset for Arbor where we have only access to the plastic
        # contribution of the weight
        if i == 4:
            weight_offset = config_lif["synapses"]["cond_exp_stdp"]["weight"]
        else:
            weight_offset = 0

        axes.flat[i - 1].plot(arbor_data[:, 0], arbor_data[:, i] + weight_offset,
                              label='Arbor', marker='None')
        axes.flat[i - 1].plot(brian_data[:, 0], brian_data[:, i],
                              label='Brian', linestyle='dashed', marker='None')
        axes.flat[i - 1].set_xlabel("Time (ms)")

    axes.flat[0].set_ylabel("Membrane potential (mV)")
    axes.flat[1].set_ylabel("Excitatory conductance (µS)")
    axes.flat[2].set_ylabel("Inhibitory conductance (µS)")
    axes.flat[3].set_ylabel("Excitatory weight (µS)")

    # plot spike comparison for LIF neurons
    arbor_spikes = np.loadtxt(f'arbor_spikes_{variant}_lif.dat')
    brian_spikes = np.loadtxt(f'brian2_spikes_{variant}_lif.dat')

    if len(arbor_spikes) == len(brian_spikes):
        axes.flat[4].plot(arbor_spikes, (arbor_spikes -
                                         brian_spikes) / brian_spikes * 100)
        axes.flat[4].set_xlabel("Time (ms)")
        axes.flat[4].set_ylabel("Spike time mismatch (%)")
    else:
        axes.flat[4].text(0.05, 0.5, "Mismatch in number of spikes:\nArbor: {}\nBrian2: {}".format(
            len(arbor_spikes), len(brian_spikes)))
        
    #axes.flat[-1].set_axis_off()
        
    # plot comparison for classical STDP curve
    arbor_classical_data = np.loadtxt(f'arbor_traces_{variant}_classical.dat')
    brian_classical_data = np.loadtxt(f'brian2_traces_{variant}_classical.dat')
    axes.flat[5].plot(arbor_classical_data[:, 0], arbor_classical_data[:, 1],
                      label='Arbor', marker='None')
    axes.flat[5].plot(brian_classical_data[:, 0], brian_classical_data[:, 1],
                      label='Brian', linestyle='dashed', marker='None')
    axes.flat[5].set_xlabel(r'Time difference $\Delta t$ (ms)')
    axes.flat[5].set_ylabel(r'Weight change $\Delta w$ (µS)')

    fig.legend(*axes[0, 0].get_legend_handles_labels(), loc="upper center")

    fig.savefig(f'comparison_{variant}.png')
    fig.savefig(f'comparison_{variant}.svg')


if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('variant', help="name of variant, e.g., brian2_arbor")
    args = parser.parse_args()

    main(args.variant)
