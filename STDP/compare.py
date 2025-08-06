#!/usr/bin/env python3
"""Plots to compare results from Arbor and Brian 2 or NEST
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
    

def main(variant):
    """
    Plots comparison of spikes and traces between Arbor and Brian 2 or NEST
    """
    config_lif = json.load(open("config_lif.json"))
    config_classical = json.load(open("config_classical.json"))

    # load data and check variant
    trace_data_lif_arbor = np.loadtxt(f'arbor_traces_{variant}_lif.dat')
    spike_data_lif_arbor = np.loadtxt(f'arbor_spikes_{variant}_lif.dat')
    trace_data_classical_arbor = np.loadtxt(f'arbor_traces_{variant}_classical.dat')
    if variant == "brian2_arbor":
        trace_data_lif_ref = np.loadtxt(f'brian2_traces_{variant}_lif.dat')
        spike_data_lif_ref = np.loadtxt(f'brian2_spikes_{variant}_lif.dat')
        trace_data_classical_ref = np.loadtxt(f'brian2_traces_{variant}_classical.dat')
        ref_name = "Brian"
    elif variant == "nest_arbor":
        trace_data_lif_ref = np.loadtxt(f'nest_traces_{variant}_lif.dat')
        spike_data_lif_ref = np.loadtxt(f'nest_spikes_{variant}_lif.dat')
        trace_data_classical_ref_without_time = np.loadtxt(f'nest_traces_{variant}_classical.dat')
        trace_data_classical_ref = np.column_stack(
          [trace_data_classical_arbor[:, 0],
           trace_data_classical_ref_without_time]) # adding times from Arbor data (no data because NEST always has one addition postsynpatic spike, see "nest_stdp_classical.py")
        ref_name = "NEST"
    else:
        raise ValueError(f"Unsupported variant: '{variant}'.")

    fig, axes = plt.subplots(nrows=3, ncols=2, sharex=False, figsize=(10, 10))

    stimulus_times_exc = config_lif["synapses"]["cond_exp_stdp"]["stimulus_times"]
    stimulus_times_inh = config_lif["synapses"]["cond_exp"]["stimulus_times"]

    print(f"Number of exc. spikes: {len(stimulus_times_exc)}")
    print(f"Number of inh. spikes: {len(stimulus_times_inh)}")

    # plot trace comparisons for LIF neurons

    print("---------------------------------\n"
          "Plotting LIF trace comparison.")

    # time, membrane voltage, excitatory conductance, inhibitory conductance,
    # weight
    ylabels = ["Membrane potential (mV)", "Excitatory conductance (µS)", "Inhibitory conductance (µS)",
               "Excitatory weight (µS)"]
    for i, ylabel in enumerate(ylabels):
        print(f"Considering {ylabel} (panel {i})")

        for spike_time in stimulus_times_exc:
            axes.flat[i].axvline(spike_time, ymin=0.9, ymax=1, color='blue')

        for spike_time in stimulus_times_inh:
            axes.flat[i].axvline(spike_time, ymin=0.85,
                                     ymax=0.95, color='red')

        # offset for Arbor where we have only access to the plastic
        # contribution of the weight
        if i == 4:
            weight_offset = config_lif["synapses"]["cond_exp_stdp"]["weight"]
        else:
            weight_offset = 0

        axes.flat[i].plot(trace_data_lif_arbor[:, 0], trace_data_lif_arbor[:, i+1] + weight_offset,
                              label='Arbor', marker='None')
        axes.flat[i].plot(trace_data_lif_ref[:, 0], trace_data_lif_ref[:, i+1],
                              label=ref_name, linestyle='dashed', marker='None')
        axes.flat[i].set_xlabel("Time (ms)")

        axes.flat[i].set_ylabel(ylabel)

        # compute and print R^2 and RMSE
        print(f"{ylabel}:")
        compute_and_print_goodness(trace_data_lif_arbor[:, i+1], trace_data_lif_ref[:, i+1])


    # plot spike comparison for LIF neurons
    if len(spike_data_lif_arbor) == len(spike_data_lif_ref):
        axes.flat[4].plot(spike_data_lif_arbor, (spike_data_lif_arbor -
                                         spike_data_lif_ref) / spike_data_lif_ref * 100)
        axes.flat[4].set_xlabel("Time (ms)")
        axes.flat[4].set_ylabel("Spike time mismatch (%)")
    else:
        axes.flat[4].text(0.05, 0.5, "Mismatch in number of spikes:\nArbor: {}\n{} {}".format(
            len(spike_data_lif_arbor), ref_name, len(spike_data_lif_ref)))
        
    #axes.flat[-1].set_axis_off()
        
    print("--------------------------------------\n"
          "Plotting classical STDP curves.")

    # plot comparison for classical STDP curve
    axes.flat[5].plot(trace_data_classical_arbor[:, 0], trace_data_classical_arbor[:, 1],
                      label='Arbor', marker='None')
    axes.flat[5].plot(trace_data_classical_ref[:, 0], trace_data_classical_ref[:, 1],
                      label=ref_name, linestyle='dashed', marker='None')
    axes.flat[5].set_xlabel(r'Time delay $\Delta t$ (ms)')
    axes.flat[5].set_ylabel(r'Weight change $\Delta w$ (µS)')

    fig.legend(*axes[0, 0].get_legend_handles_labels(), loc="lower right")

    fig.tight_layout()
    fig.savefig(f'comparison_{variant}.png')
    fig.savefig(f'comparison_{variant}.svg')

    # compute and print R^2 and RMSE
    compute_and_print_goodness(trace_data_classical_arbor[:, 1], trace_data_classical_ref[:, 1])



if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('variant', help="name of variant, e.g., brian2_arbor")
    args = parser.parse_args()

    main(args.variant)
