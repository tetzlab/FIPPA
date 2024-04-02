#!/usr/bin/env python3
"""Brian 2 simulation of a single cell receiving inhibitory and plastic
excitatory stimulus.
"""


import json
import numpy as np
from brian2 import siemens, cmeter, farad, mV, meter, umeter, ms, uS
from brian2 import SpikeGeneratorGroup, NeuronGroup, Synapses, StateMonitor, SpikeMonitor
from brian2 import run, defaultclock


def main(variant):
    """Runs simulation with spikes generated at specific times
    and stores results."""

    config = json.load(open(f"config_{variant}_lif.json"))
    neuron_config = config["neuron"]
    # sphere with 200 um radius
    area = 4 * np.pi * (neuron_config["radius"] * umeter)**2

    C_mem = area * \
        neuron_config["specific_capacitance"] * farad / (meter * meter)
    gl = area * neuron_config["g_leak"] * siemens / (cmeter * cmeter)

    syn_config_stdp = config["synapses"]["cond_exp_stdp"]
    syn_config = config["synapses"]["cond_exp"]

    e_leak = neuron_config["e_leak"] * mV
    e_exc = syn_config_stdp["reversal_potential"] * mV
    e_inh = syn_config["reversal_potential"] * mV
    v_thresh = neuron_config["v_thresh"] * mV
    v_reset = neuron_config["v_reset"] * mV
    tau_e = syn_config_stdp["tau"] * ms
    tau_i = syn_config["tau"] * ms
    tau_refrac = neuron_config["tau_refrac"] * ms

    defaultclock.dt = config["simulation"]["dt"] * ms

    eqs_neurons = '''
    dv/dt = gl/C_mem * (e_leak - v) + ge/C_mem*(e_exc - v) + gi/C_mem*(e_inh-v): volt (unless refractory)
    dge/dt = -ge / tau_e : siemens
    dgi/dt = -gi / tau_i : siemens
    '''

    spike_times_exc = np.array(syn_config_stdp["stimulus_times"]) * ms
    spike_times_inh = np.array(syn_config["stimulus_times"]) * ms

    ssg_exc = SpikeGeneratorGroup(
        1, [0] * len(spike_times_exc), spike_times_exc)
    ssg_inh = SpikeGeneratorGroup(
        1, [0] * len(spike_times_inh), spike_times_inh)

    neurons = NeuronGroup(1, eqs_neurons, threshold='v > v_thresh',
                          reset='v = v_reset', refractory=tau_refrac, method='euler')
    neurons.v = e_leak

    tau_pre = syn_config_stdp["tau_pre"] * ms
    tau_post = syn_config_stdp["tau_post"] * ms
    A_pre = syn_config_stdp["A_pre"] * uS
    A_post = syn_config_stdp["A_post"] * uS

    S_exc = Synapses(ssg_exc, neurons,
                 '''w : siemens
                 dapre/dt = -apre/tau_pre : siemens (event-driven)
                 dapost/dt = -apost/tau_post : siemens (event-driven)
                 ''',
                 on_pre='''
                 ge += w
                 apre += A_pre
                 w += apost
                 ''',
                 on_post='''
                 apost += A_post
                 w += apre
                 ''',
                 delay=0 * ms)

    S_exc.connect('True')
    S_exc.w = syn_config_stdp["weight"] * uS

    S_inh = Synapses(ssg_inh, neurons, '''w : siemens''',
                     on_pre='''gi += w''', delay=0 * ms)
    S_inh.connect('True')
    S_inh.w = syn_config["weight"] * uS

    neuron_monitor = StateMonitor(neurons, ['v', 'ge', 'gi'], record=True)
    synapse_monitor = StateMonitor(S_exc, ['w'], record=True)
    spikemon = SpikeMonitor(neurons)

    run(config["simulation"]["runtime"] * ms)

    np.savetxt(f'brian2_traces_{variant}_lif.dat', np.column_stack(
        [neuron_monitor.t / ms, neuron_monitor.v[0] / mV,
         neuron_monitor.ge[0] / uS, neuron_monitor.gi[0] / uS,
         synapse_monitor.w[0] / uS]))
    np.savetxt(f'brian2_spikes_{variant}_lif.dat', spikemon.t / ms)


if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('variant')
    args = parser.parse_args()

    main(args.variant)
