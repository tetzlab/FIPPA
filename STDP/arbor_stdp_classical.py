#!/usr/bin/env python3
"""
Arbor simulation of two neuron populations connecting via STDP synapses.

Event generators are not used; instead, the spiking is inherently triggered in mechanisms,
resembling the way of the Brian 2 implementation in 'brian2_stdp_classical.py'.
"""

import json
import arbor
import numpy as np


class SingleRecipe(arbor.recipe):
    """Implementation of Arbor simulation recipe."""

    def __init__(self, config):
        """Initialize the recipe from config."""

        # The base C++ class constructor must be called first, to ensure that
        # all memory in the C++ class is initialized correctly.
        arbor.recipe.__init__(self)

        self.the_props = arbor.neuron_cable_properties()
        self.the_cat = arbor.load_catalogue("./custom-catalogue.so")
        self.the_cat.extend(arbor.default_catalogue(), "")
        self.the_props.catalogue = self.the_cat

        self.config = config
        self.N = config["simulation"]["N"]
        self.dt = self.config["simulation"]["dt"]
        self.t_max = self.config["simulation"]["runtime"]

        # arrays of spike time values
        self.t_spike_1 = np.array([ ])
        self.t_spike_2 = np.array([ ])


    def num_cells(self):
        """Return the number of cells."""
        return 2*self.N


    def num_sources(self, gid):
        """Return the number of spikes sources on gid."""
        if gid < self.N:
            return 0
        else:
            return 1


    def num_targets(self, gid):
        """Return the number of post-synaptic targets on gid."""
        if gid < self.N:
            return 1
        else:
            return 0


    def cell_kind(self, gid):
        """Return type of cell with gid."""
        return arbor.cell_kind.cable


    def cell_description(self, gid):
        """Return cell description of gid."""

        # morphology
        tree = arbor.segment_tree()
        radius = self.config["neuron"]["radius"]

        tree.append(arbor.mnpos,
                    arbor.mpoint(-radius, 0, 0, radius),
                    arbor.mpoint(radius, 0, 0, radius),
                    tag=1)

        labels = arbor.label_dict({'center': '(location 0 0.5)'})

        # cell mechanism
        e_thresh = self.the_cat[self.config["neuron"]["type"]].parameters["e_thresh"].default
        e_reset = self.the_cat[self.config["neuron"]["type"]].parameters["e_reset"].default
        decor = arbor.decor()
        decor.set_property(Vm=e_reset)
        neuron = arbor.mechanism(self.config["neuron"]["type"])
        neuron.set("tau_refrac", self.config["neuron"]["tau_refrac"])
        if gid < self.N:
            # define spike times for neurons 0 to N-1
            t_spike = gid*self.t_max/(self.N-1)
            neuron.set("t_spike", t_spike)
            try:
                self.t_spike_1 = np.column_stack((self.t_spike_1, [t_spike, gid]))
            except ValueError:
                self.t_spike_1 = [t_spike, gid]
        else:
            # define spike times for neurons N to 2*N
            t_spike = (2*self.N-1-gid)*self.t_max/(self.N-1)
            neuron.set("t_spike", t_spike)
            try:
                self.t_spike_2 = np.column_stack((self.t_spike_2, [t_spike, gid]))
            except ValueError:
                self.t_spike_2 = [t_spike, gid]

            # add incoming plastic synapse
            syn_config_stdp = self.config["synapses"]["stdp"]
            mech_expsyn = arbor.mechanism('expsyn_stdp')
            mech_expsyn.set('taupre', syn_config_stdp["tau_pre"])
            mech_expsyn.set('taupost', syn_config_stdp["tau_post"])
            mech_expsyn.set('Apre', syn_config_stdp["A_pre"])
            mech_expsyn.set('Apost', syn_config_stdp["A_post"])
            mech_expsyn.set('max_weight', 50)
            decor.place('"center"', arbor.synapse(mech_expsyn), "expsyn_stdp_post")

        decor.place('"center"', arbor.threshold_detector(e_thresh), "spike_detector")
        decor.paint('(all)', arbor.density(neuron))

        return arbor.cable_cell(tree, decor, labels)


    def connections_on(self, gid):
        """Defines the list of synaptic connections incoming to the neuron given by gid"""
        
        policy = arbor.selection_policy.univalent
        weight = 0
        delay = self.dt # may not be <= 0

        # neurons with gid 0 to N-1 are presynaptic
        if gid < self.N:
            conn = [ ]
        
        # neurons with gid N to 2*N are postsynaptic
        else:
            src = gid - self.N
            conn = [arbor.connection((src, "spike_detector"), ('expsyn_stdp_post', policy), weight, delay)]

        return conn


    def probes(self, gid):
        """Return probes on gid."""

        probe_list = []
        #probe_list = [arbor.cable_probe_membrane_voltage('"center"')]
        #probe_list = [arbor.cable_probe_density_state('"center"', self.config["neuron"]["type"], "t")]

        # neurons with gid N to 2*N are postsynaptic
        if gid >= self.N and gid < 2*self.N:
            probe_list.append(arbor.cable_probe_point_state(0, "expsyn_stdp", "weight_plastic"))
        
        return probe_list


    def global_properties(self, kind):
        """Return the global properties."""
        assert kind == arbor.cell_kind.cable

        return self.the_props


def main(variant):
    """Runs simulation and stores results."""

    # set up simulation and run
    config = json.load(open("config_classical.json"))
    recipe = SingleRecipe(config)

    context = arbor.context()
    domains = arbor.partition_load_balance(recipe, context)
    sim = arbor.simulation(recipe, context, domains)

    sim.record(arbor.spike_recording.all)
    reg_sched = arbor.regular_schedule(config["simulation"]["dt"])
    handle_weight_plastic_array = [sim.sample((i, 0), reg_sched) for i in range(recipe.N, 2*recipe.N)]

    sim.run(tfinal=config["simulation"]["runtime"] + 1,
            dt=config["simulation"]["dt"])

    # read out and store weight changes and spike data
    data_weight_plastic_final = np.zeros(recipe.N)
    for i in range(recipe.N):
        if len(sim.samples(handle_weight_plastic_array[i])) > 0:
            data_buf, _ = sim.samples(handle_weight_plastic_array[i])[0]
            data_weight_plastic_final[i] = data_buf[-1, 1]

    t_spike_1_unsorted_T = recipe.t_spike_1.T
    t_spike_2_unsorted_T = recipe.t_spike_2.T
    t_spike_1 = t_spike_1_unsorted_T[t_spike_1_unsorted_T[:,1].argsort()].T
    t_spike_2 = t_spike_2_unsorted_T[t_spike_2_unsorted_T[:,1].argsort()].T

    data_stacked = np.column_stack(
        [t_spike_2[0] - t_spike_1[0],
         data_weight_plastic_final])
    
    spikes = np.column_stack((sim.spikes()['time'], sim.spikes()['source']['gid']))

    np.savetxt(f'arbor_traces_{variant}_classical.dat', data_stacked)
    np.savetxt(f'arbor_spikes_{variant}_classical.dat', spikes, fmt="%.4f %.0f") # integer formatting for neuron number
    

if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('variant', help="name of variant, e.g., brian2_arbor")
    args = parser.parse_args()

    main(args.variant)
