#!/usr/bin/env python3
"""Brian 2 simulation of a single cell receiving inhibitory and plastic
excitatory stimulus (from https://brian2.readthedocs.io/en/stable/examples/synapses.spike_based_homeostasis.html).
"""

import itertools
import json
import numpy as np
import matplotlib.pyplot as plt

from brian2 import TimedArray, PoissonGroup, NeuronGroup, Synapses
from brian2 import StateMonitor, SpikeMonitor, PopulationRateMonitor
from brian2 import defaultclock, run
from brian2 import Hz, ms, second, mV, nA, mS

def main(config_file, considered_neuron = 1):
    """Runs simulation and stores results.
    
    Parameters
    ----------
    config_file
        Name of the JSON configuration file.
    considered_neuron
        The neuron to monitor: 0 for neuron without homeostasis, 1 for neuron with
        homeostasis.
    """

    config = json.load(open(config_file, 'r'))

    # The synaptic weight from the steady stimulus is plastic
    syn_config_steady = config["stimulus"]["steady"]
    steady_stimulus = TimedArray(syn_config_steady["rates"]*Hz, dt=syn_config_steady["dt"]*ms)
    steady_poisson = PoissonGroup(1, rates='steady_stimulus(t)')

    # The synaptic weight from the varying stimulus is static
    syn_config_varying = config["stimulus"]["varying"]
    varying_stimulus = TimedArray(syn_config_varying["rates"]*Hz, dt=syn_config_varying["dt"]*ms)
    varying_poisson = PoissonGroup(1, rates='varying_stimulus(t)')

    # dw_plus/dw_minus determines scales the steady stimulus rate to the target firing rate, must not be larger 1
    # the magntude of dw_plus and dw_minus determines the "speed" of the homeostasis
    parameters = {
        'tau': config["neuron"]["tau_e"]*ms,  # membrane time constant
        't_ref': config["neuron"]["t_ref"]*ms,  # neuronal refractory period
        'v_thresh': config["neuron"]["v_thresh"]*mV,  # neuronal threshold potential
        'v_reset': config["neuron"]["e_reset"]*mV, # neuronal reset potential
        'v_rev': config["neuron"]["e_leak"]*mV, # neuronal reversal potential
        'g_leak': config["neuron"]["g_leak"]*mS,  # membrane conductance
        'dw_plus': syn_config_steady["dw_plus"]*nA,  # weight increment on presynaptic spike
        'dw_minus': syn_config_steady["dw_minus"]*nA,  # weight increment on postsynaptic spike
        'w_init': syn_config_steady["w_init"]*nA,  # initial plastic weight
        'w_max': syn_config_steady["w_max"]*nA,  # maximum plastic weight
        'w_non_plastic': syn_config_varying["w_non_plastic"]*nA  # weight of the non-plastic synapse
    }

    eqs = 'dv/dt = (v_rev - v)/tau : volt (unless refractory)'
    g_leak = parameters['g_leak']

    # Target neuron with homeostasis
    neuron_with_homeostasis = NeuronGroup(1, eqs,
                                         threshold='v > v_thresh', reset='v = v_reset',
                                         method='euler', refractory='t_ref',
                                         namespace=parameters)
    neuron_with_homeostasis.v = parameters['v_rev']

    # Target neuron without homeostasis
    neuron_without_homeostasis = NeuronGroup(1, eqs,
                                             threshold='v > v_thresh', reset='v = v_reset',
                                             method='euler', refractory='t_ref',
                                             namespace=parameters)
    neuron_without_homeostasis.v = parameters['v_rev']

    # Steady-rate input synapse (plastic)
    plastic_synapse = Synapses(steady_poisson, neuron_with_homeostasis,
                               'w : amp',
                               on_pre='''
                               v_post += w / g_leak
                               w = clip(w + dw_plus, 0*namp, w_max)
                               ''',
                               on_post='''
                               w = clip(w - dw_minus, 0*namp, w_max)
                               ''', namespace=parameters)
    plastic_synapse.connect()
    plastic_synapse.w = parameters['w_init']

    # Steady-rate input synapse (non-plastic)
    non_plastic_synapse_neuron_without_homeostasis = Synapses(varying_poisson,
                                                              neuron_without_homeostasis,
                                                              'w : amp', on_pre='v_post += w / g_leak')
    non_plastic_synapse_neuron_without_homeostasis.connect()
    non_plastic_synapse_neuron_without_homeostasis.w = parameters['w_non_plastic']

    non_plastic_synapse_neuron_with_homeostasis = Synapses(varying_poisson,
                                                           neuron_with_homeostasis,
                                                           'w : amp', on_pre='v_post += w / g_leak')
    non_plastic_synapse_neuron_with_homeostasis.connect()
    non_plastic_synapse_neuron_with_homeostasis.w = parameters['w_non_plastic']

    M_v_with = StateMonitor(neuron_with_homeostasis, 'v', record=True)
    M_v_without = StateMonitor(neuron_without_homeostasis, 'v', record=True)
    M_w = StateMonitor(plastic_synapse, 'w', record=True)
    #M_rate_neuron_with_homeostasis = PopulationRateMonitor(neuron_with_homeostasis)
    #M_rate_neuron_without_homeostasis = PopulationRateMonitor(neuron_without_homeostasis)
    SM_with = SpikeMonitor(neuron_with_homeostasis)
    SM_without = SpikeMonitor(neuron_without_homeostasis)

    duration = config["simulation"]["runtime"]*ms
    defaultclock.dt = config["simulation"]["dt"]*ms
    run(duration, report='text')

    # collect and store data
    input_dts = np.arange(0., len(varying_stimulus.values)*varying_stimulus.dt/ms, varying_stimulus.dt/ms)
    input_x = list(itertools.chain(*zip(input_dts, input_dts)))
    input_y = [0] + list(itertools.chain(*zip(varying_stimulus.values, varying_stimulus.values)))[:-1]
    input_stacked = np.column_stack([input_x, input_y])
    if considered_neuron == 1:
        results_stacked = np.column_stack([M_v_with.t/ms, M_v_with.v[0]/mV, M_w.w[0]/nA])
        spike_times = SM_with.t/ms
    else:
        results_stacked = np.column_stack([M_v_without.t/ms, M_v_without.v[0]/mV, np.zeros_like(M_v_with.v[0])])
        spike_times = SM_without.t/ms
    np.savetxt(f'brian2_input_{considered_neuron}.dat', input_stacked)

    return results_stacked, spike_times


if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('config', help="name of config file")
    parser.add_argument('num_trials', type=int, help="number of trials to consider")
    parser.add_argument('considered_neuron', type=int, help="the type of neuron to consider")
    args = parser.parse_args()

    num_trials = args.num_trials
    data_stacked_sum = np.array([])
    spike_times_all = []
    for i in range(num_trials):
        data_stacked, spike_times = main(args.config, args.considered_neuron)
        if data_stacked_sum.size == 0:
            data_stacked_sum = np.array(data_stacked)
        else:
            data_stacked_sum += data_stacked
        spike_times_all.extend(spike_times)
    
    np.savetxt(f'brian2_traces_{args.considered_neuron}.dat', data_stacked_sum / num_trials)
    np.savetxt(f'brian2_spikes_{args.considered_neuron}.dat', spike_times_all)
