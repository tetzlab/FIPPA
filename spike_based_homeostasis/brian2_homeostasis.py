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
from brian2 import Hz, ms, second, mV

def main(config_file):
    """Runs simulation and stores results."""

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
        'v_reset': config["neuron"]["v_reset"]*mV, # neuronal reset potential
        'dw_plus': syn_config_steady["dw_plus"]*mV,  # weight increment on pre spike
        'dw_minus': syn_config_steady["dw_minus"]*mV,  # weight increment on post spike
        'w_initial': 0*mV,  # initial plastic weight
        'w_max': syn_config_steady["w_max"]*mV,  # maximum plastic weight
        'w_non_plastic': syn_config_varying["weight"]*mV  # weight of the non-plastic synapse
    }

    eqs = 'dv/dt = (- v)/tau : volt (unless refractory)'

    neuron_with_homeostasis = NeuronGroup(1, eqs,
                                        threshold='v > v_thresh', reset='v = v_reset',
                                        method='euler', refractory='t_ref',
                                        namespace=parameters)
    neuron_without_homeostasis = NeuronGroup(1, eqs,
                                            threshold='v > v_thresh', reset='v = v_reset',
                                            method='euler', refractory='t_ref',
                                            namespace=parameters)

    plastic_synapse = Synapses(steady_poisson, neuron_with_homeostasis,
                            'w : volt',
                            on_pre='''
                            v_post += w
                            w = clip(w + dw_plus, 0*mvolt, w_max)
                            ''',
                            on_post='''
                            w = clip(w - dw_minus, 0*mvolt, w_max)
                            ''', namespace=parameters)
    plastic_synapse.connect()
    plastic_synapse.w = parameters['w_initial']

    non_plastic_synapse_neuron_without_homeostasis = Synapses(varying_poisson,
                                                            neuron_without_homeostasis,
                                                            'w : volt', on_pre='v_post += w')
    non_plastic_synapse_neuron_without_homeostasis.connect()
    non_plastic_synapse_neuron_without_homeostasis.w = parameters['w_non_plastic']

    non_plastic_synapse_neuron = Synapses(varying_poisson, neuron_with_homeostasis,
                                        'w : volt', on_pre='v_post += w')
    non_plastic_synapse_neuron.connect()
    non_plastic_synapse_neuron.w = 2*mV

    M_v_with = StateMonitor(neuron_with_homeostasis, 'v', record=True)
    #M_v_without = StateMonitor(neuron_without_homeostasis, 'v', record=True)
    M_w = StateMonitor(plastic_synapse, 'w', record=True)
    #M_rate_neuron_with_homeostasis = PopulationRateMonitor(neuron_with_homeostasis)
    #M_rate_neuron_without_homeostasis = PopulationRateMonitor(neuron_without_homeostasis)
    SM_with = SpikeMonitor(neuron_with_homeostasis)
    SM_without = SpikeMonitor(neuron_without_homeostasis)

    duration = config["simulation"]["runtime"]*ms
    defaultclock.dt = config["simulation"]["dt"]*ms
    run(duration, report='text')

    # collect data and store
    input_dts = np.arange(0., len(varying_stimulus.values)*varying_stimulus.dt/ms, varying_stimulus.dt/ms)
    input_x = list(itertools.chain(*zip(input_dts, input_dts)))
    input_y = [0] + list(itertools.chain(*zip(varying_stimulus.values, varying_stimulus.values)))[:-1]
    input_stacked = np.column_stack([input_x, input_y])
    results_stacked = np.column_stack([M_v_with.t/ms, M_v_with.v[0]/mV, M_w.w[0]/mV])
    spike_times_with = SM_with.t/ms
    spike_times_without = SM_without.t/ms
    np.savetxt(f'brian2_input.dat', input_stacked)
    np.savetxt(f'brian2_traces.dat', results_stacked)
    np.savetxt(f'brian2_spikes.dat', spike_times_with)
    np.savetxt(f'brian2_spikes_without_homestasis.dat', spike_times_without)


if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('config', help="name of config file")
    args = parser.parse_args()

    main(args.config)