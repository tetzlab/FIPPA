#!/usr/bin/env python3
"""NEST simulation of a single cell receiving inhibitory and plastic
excitatory stimulus.

TODO The implementation is not complete yet, some things need to be adapted or added.
"""

import nest
import json
import numpy as np

def main(variant):
    """Runs simulation with spikes generated at specific times
    and stores results."""

    # set parameters
    config = json.load(open("config_lif.json"))
    neuron_config = config["neuron"]
    area = 4 * np.pi * (neuron_config["radius"])**2 # in um^2
    C_mem = area * \
        neuron_config["specific_capacitance"]
    gl = area * neuron_config["g_leak"]
    syn_config_stdp = config["synapses"]["cond_exp_stdp"] # config for excitatory, plastic synapse
    syn_config = config["synapses"]["cond_exp"] # config for inhibitory, static synapse
    e_leak = neuron_config["e_leak"]
    e_exc = syn_config_stdp["reversal_potential"]
    e_inh = syn_config["reversal_potential"]
    v_thresh = neuron_config["v_thresh"]
    v_reset = neuron_config["v_reset"]
    tau_exc = syn_config_stdp["tau"]
    tau_inh = syn_config["tau"]
    t_ref = neuron_config["tau_refrac"]
    tau_pre = syn_config_stdp["tau_pre"]
    tau_post = syn_config_stdp["tau_post"]
    A_pre = syn_config_stdp["A_pre"]
    A_post = syn_config_stdp["A_post"]
    tmax = config["simulation"]["runtime"]  # simulation time in ms
    res = config["simulation"]["dt"] # temporal resolution in ms

    # reset NEST kernel and set timestep resolution
    nest.ResetKernel()
    nest.SetKernelStatus({"resolution": res})
    
    # define neuron parameters
    neuron_params = {
        'V_m': e_leak,    # Initial membrane ptoential in mV
        'C_m': C_mem,      # Membrane capacitance in pF
        'g_L': gl,     # Leak conductance in nS
        'E_L': e_leak,    # Leak reversal potential in mV
        'V_th': v_thresh,   # Threshold potential in mV
        'V_reset': v_reset,# Reset potential in mV
        't_ref': t_ref,    # Absolute refractory period in ms
        'tau_syn_ex': tau_exc,   # Exponential decay time constant of excitatory synaptic conductance kernel in ms
        'tau_syn_in': tau_inh,   # Exponential decay time constant of inhibitory synaptic conductance kernel in ms
        'E_ex': e_exc,        # Excitatory reversal potential in mV
        'E_in': e_inh,        # Inhibitory reversal potential in mV
        'I_e': 0.0,       # External current in pA
        'tau_minus' : tau_post    # Time constant for STDP depression in ms
    }

    # create neurons
    relay_neuron = nest.Create('parrot_neuron', 1)
    neuron = nest.Create('iaf_cond_exp', 1, params=neuron_params)

    # load spike times from configuration and create spike generators
    spike_times_exc = np.array(syn_config_stdp["stimulus_times"])  # in ms
    spike_times_inh = np.array(syn_config["stimulus_times"])      # in ms
    sg_exc = nest.Create('spike_generator', 1, {'spike_times': np.array(spike_times_exc, dtype=float)})
    sg_inh = nest.Create('spike_generator', 1, {'spike_times': np.array(spike_times_inh, dtype=float)})

    # define STDP synapse parameters (cf. https://nest-simulator.readthedocs.io/en/v3.8/models/stdp_synapse.html)
    stdp_params = {
        'synapse_model': 'stdp_synapse',
        'weight': syn_config_stdp["weight"],                 # initial weight in uS
        'delay': res,                  # synaptic delay in ms
        'lambda': res,                 # step size/learning rate in ms
        'alpha': np.abs(A_post)/A_pre,   # asymmetry parameter (A_{-} = alpha*A_{+})
        'mu_plus': 0.0,                # zero for additive STDP
        'mu_minus': 0.0,               # zero for additive STDP
        'tau_plus': tau_pre,            # Time constant for potentiation (tau_minus is defined in postsynaptic neuron) in ms
        'Wmax' : np.max((np.abs(A_post), A_pre)) + 1.0    # maximum weight
    }

    # create excitatory (plastic) and inhibitory (static) synapse; the plastic synapse requires a relay neuron
    nest.Connect(sg_exc, relay_neuron)
    nest.Connect(relay_neuron,  neuron, syn_spec=stdp_params)
    nest.Connect(sg_inh, neuron, syn_spec={"weight" : syn_config["weight"]})

    # create multimeter to record membrane potential, excitatory and inhibitory conductances
    multimeter = nest.Create('multimeter', 1, {
        'interval': 0.01,  # record every 0.01 ms
        'record_from': ['V_m', 'g_ex', 'g_in']
    })

    # create spike recorder
    spike_recorder = nest.Create('spike_recorder')

    # Connect monitors
    nest.Connect(multimeter, neuron)
    nest.Connect(neuron, spike_recorder)

    # Run simulation
    nest.Simulate(tmax+res)

    # Get traces data and save them
    traces_data = nest.GetStatus(multimeter, 'events')[0]
    np.savetxt(f'nest_traces_{variant}_lif.dat',
            np.column_stack([traces_data['times'],
                             traces_data['V_m'],
                             traces_data['g_ex'],
                             traces_data['g_in'],
                             traces_data['g_in']]))

    # Get spike data and save them
    spike_data = nest.GetStatus(spike_recorder, 'events')[0]
    np.savetxt(f'nest_spikes_{variant}_lif.dat',  
               np.sort(np.column_stack([spike_data['times'], 
                                        spike_data['senders']]), axis=0),
               fmt="%.4f %.0f") # integer formatting for neuron number

if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('variant')
    args = parser.parse_args()

    main(args.variant)
