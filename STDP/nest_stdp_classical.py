#!/usr/bin/env python3
"""
NEST simulation of two neuron populations connecting via STDP synapses.
Generates classical STDP curve.
"""
import nest
import json
import numpy as np

def main(variant):
    """Runs simulation of classical STDP curve (implementation inspired by 
    https://brian2.readthedocs.io/en/stable/resources/tutorials/2-intro-to-brian-synapses.html)
    and stores results."""

    # set parameters
    config = json.load(open("config_classical.json"))
    taupre = config["synapses"]["stdp"]["tau_pre"] # ms
    taupost = config["synapses"]["stdp"]["tau_post"] # ms
    Apre = config["synapses"]["stdp"]["A_pre"] # uS
    Apost = -Apre*taupre/taupost*1.05 # uS
    tmax = config["simulation"]["runtime"] # runtime of STDP window simulation in ms
    tmaxall = 5*tmax # total runtime in ms (since potentiation requires a second, dummy presynaptic spike)
    N = config["simulation"]["N"]
    t_ref = config["neuron"]["tau_refrac"] # ms
    res = config["simulation"]["dt"] # temporal resolution in ms

    # reset NEST kernel and set timestep resolution
    nest.ResetKernel()
    nest.SetKernelStatus({"resolution": res})

    # create neurons (pre just relays spikes, post is full 'iaf_psc_alpha' neuron with refractory period such that it spikes only once)
    pre_neuron = nest.Create('parrot_neuron', N)
    post_neuron = nest.Create('iaf_psc_alpha', N, params={"t_ref": t_ref, "tau_syn_ex" : res, 'tau_minus' : taupost})

    # create spike generators for pre- and postsynaptic neurons
    spike_times_pre = np.around([i*tmax/(N-1)+res for i in range(N)], int(-np.log10(res))) # rounding because decimal places need to match the resolution set above
    spike_times_post = np.around([(N-1-i)*tmax/(N-1)+res for i in range(N)], int(-np.log10(res))) # rounding because decimal places need to match the resolution set above
    spike_generator_pre = nest.Create('spike_generator', N) #, params={'spike_times': spike_times_pre})
    spike_generator_post = nest.Create('spike_generator', N) #, params={'spike_times': spike_times_post})

    # go through the pre- and postsynaptic generators and set the desired spike times
    for i in range(N):
        # for depression - one presynaptic spike suffices
        if i >= N // 2:
            nest.SetStatus(nest.NodeCollection(spike_generator_pre[i]), {"spike_times": [spike_times_pre[i]]})
        # for potentiation - at least two presynaptic spikes are necessary to induce any potentiation 
        #   (cf. Supplementary of https://doi.org/10.3389/fninf.2024.1446620)
        else:
            nest.SetStatus(nest.NodeCollection(spike_generator_pre[i]), {"spike_times": [spike_times_pre[i], 4*tmax]})
        # always use one postsynaptic spike
        nest.SetStatus(nest.NodeCollection(spike_generator_post[i]), {"spike_times": [spike_times_post[i]]})

    # connect spike generators for pre- and postsynaptic neurons (unless parrot neuron, choose very high weight to enable definite spiking)
    nest.Connect(spike_generator_pre, pre_neuron, 'one_to_one')
    nest.Connect(spike_generator_post, post_neuron, 'one_to_one', syn_spec={"weight" : 2e5})

    # create spike recorders
    spike_recorder_pre = nest.Create('spike_recorder')
    spike_recorder_post = nest.Create('spike_recorder')

    # connect neurons with STDP synapses (cf. https://nest-simulator.readthedocs.io/en/v3.8/models/stdp_synapse.html)
    stdp_synapse = {
        'synapse_model': 'stdp_synapse',
        'weight': 1.0,                 # initial weight (set here to one since NEST cannot deal with zero weight - is later subtracted)
        'delay': res,                  # synaptic delay
        'lambda': res,                 # step size/learning rate
        'alpha': np.abs(Apost)/Apre,   # asymmetry parameter (A_{-} = alpha*A_{+})
        'mu_plus': 0.0,                # zero for additive STDP
        'mu_minus': 0.0,               # zero for additive STDP
        'tau_plus': taupre,            # STDP time constant (tau_minus is defined in postsynaptic neuron)
        'Wmax' : np.max((np.abs(Apost), Apre)) + 1.0    # maximum weight
    }
    nest.Connect(pre_neuron, post_neuron, 'one_to_one', syn_spec=stdp_synapse)

    # connect pre- and postsynaptic neurons to spike recorders
    nest.Connect(pre_neuron, spike_recorder_pre)
    nest.Connect(post_neuron, spike_recorder_post)

    # simulate
    nest.Simulate(tmaxall+1)

    # get final weights
    connections = nest.GetConnections(pre_neuron, post_neuron)
    final_weights = nest.GetStatus(connections, keys='weight')
    print(f"Final synaptic weights (total number: {len(final_weights)}): {final_weights}")

    # get spikes
    spike_data_pre = nest.GetStatus(spike_recorder_pre, 'events')[0]
    spike_data_post = nest.GetStatus(spike_recorder_post, 'events')[0]

    #print(f"Number of presynaptic spikes: {len(spike_data_pre['times'])}")
    #print(f"Number of postsynaptic spikes: {len(spike_data_post['times'])}")

    # store traces and spikes
    np.savetxt(f'nest_traces_{variant}_classical.dat',
               np.column_stack([(np.array(final_weights) - 1.0)]))
    np.savetxt(f'nest_spikes_{variant}_classical.dat', 
               np.sort(np.column_stack([np.concatenate((spike_data_pre['times'], spike_data_post['times'])), 
                                        np.concatenate((spike_data_pre['senders'], spike_data_post['senders']))]), axis=0),
               fmt="%.4f %.0f") # integer formatting for neuron number

if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('variant')
    args = parser.parse_args()

    main(args.variant)
