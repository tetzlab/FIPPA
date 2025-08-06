#!/usr/bin/env python3
"""
Brian 2 simulation of two neuron populations connecting via STDP synapses.
"""

import json
import numpy as np
from brian2 import ms, siemens, uS
from brian2 import NeuronGroup, Synapses, SpikeMonitor
from brian2 import run, defaultclock


def main(variant):
    """Runs simulation of classical STDP curve (based on 
    https://brian2.readthedocs.io/en/stable/resources/tutorials/2-intro-to-brian-synapses.html)
    and stores results."""

    config = json.load(open("config_classical.json"))
    neuron_config = config["neuron"]
    #start_scope()

    syn_config_stdp = config["synapses"]["stdp"]

    tau_refrac = neuron_config["tau_refrac"] * ms

    defaultclock.dt = config["simulation"]["dt"] * ms

    tau_pre = syn_config_stdp["tau_pre"] * ms
    tau_post = syn_config_stdp["tau_post"] * ms
    A_pre = syn_config_stdp["A_pre"] * uS
    A_post = - 1.05 * A_pre * tau_pre / tau_post
    t_max = config["simulation"]["runtime"]*ms
    N = config["simulation"]["N"]

    # Presynaptic neurons (`neurons_1`) spike at times from 0 to t_max
    # Postsynaptic neurons (`neurons_2`) spike at times from t_max to 0
    # So difference in spike times will vary from -t_max to +t_max
    neurons_1 = NeuronGroup(N, 't_spike : second', threshold='t > t_spike', refractory=tau_refrac)
    neurons_2 = NeuronGroup(N, 't_spike : second', threshold='t > t_spike', refractory=tau_refrac)
    neurons_1.t_spike = 'i*t_max/(N-1)'
    neurons_2.t_spike = '(N-1-i)*t_max/(N-1)'

    S = Synapses(neurons_1, neurons_2,
                '''
                w : siemens
                dapre/dt = -apre/tau_pre : siemens (event-driven)
                dapost/dt = -apost/tau_post : siemens (event-driven)
                ''',
                on_pre='''
                apre += A_pre
                w = w+apost
                ''',
                on_post='''
                apost += A_post
                w = w+apre
                ''')
    S.connect(j='i') # as many synapses as neurons in each group
    S.w = syn_config_stdp["weight"] * uS

    spikemon_1 = SpikeMonitor(neurons_1)
    spikemon_2 = SpikeMonitor(neurons_2)

    run(t_max + 1 * ms)

    np.savetxt(f'brian2_traces_{variant}_classical.dat',
               np.column_stack([(neurons_2.t_spike - neurons_1.t_spike) / ms, S.w / uS]))
    
    spike_indices = np.vstack((spikemon_1.i, spikemon_2.i)).flatten()
    spike_times = np.vstack((spikemon_1.t / ms, spikemon_2.t / ms)).flatten()
    spike_data = np.column_stack([spike_times, spike_indices])
    sorted_indices = np.lexsort((spike_data[:, 1], spike_data[:, 0]))
    sorted_spike_data = spike_data[sorted_indices]
    np.savetxt(f'brian2_spikes_{variant}_classical.dat', 
               sorted_spike_data,
               fmt="%.4f %.0f") # integer formatting for neuron number
    

if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('variant')
    args = parser.parse_args()

    main(args.variant)
