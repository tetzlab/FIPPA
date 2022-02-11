#!/usr/bin/env python3
"""Arbor simulation of a single cell

Arbor simulation of a single cell receiving inhibitory and plastic
excitatory stimulus.
"""

import json
import arbor
import numpy


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
        self.the_props.register(self.the_cat)

        self.config = config

    def num_cells(self):
        """Return the number of cells."""
        return 1

    def num_sources(self, gid):
        """Return the number of spikes sources on gid."""
        assert gid == 0
        return 1

    def num_targets(self, gid):
        """Return the number of post-synaptic targets on gid."""
        assert gid == 0
        return 2

    def cell_kind(self, gid):
        """Return type of cell with gid."""
        assert gid == 0
        return arbor.cell_kind.cable

    def cell_description(self, gid):
        """Return cell description of gid."""
        assert gid == 0

        neuron_config = self.config["neuron"]

        # morphology
        tree = arbor.segment_tree()
        radius = neuron_config["radius"]

        tree.append(arbor.mnpos,
                    arbor.mpoint(-radius, 0, 0, radius),
                    arbor.mpoint(radius, 0, 0, radius),
                    tag=1)

        labels = arbor.label_dict({'center': '(location 0 0.5)'})

        # cell mechanism
        decor = arbor.decor()
        decor.set_property(Vm=neuron_config["e_leak"], cm=neuron_config["specific_capacitance"])
        lif = arbor.mechanism(neuron_config["type"])
        v_thresh = neuron_config["v_thresh"]
        lif.set("e_thresh", v_thresh)
        lif.set("e_reset", neuron_config["v_reset"])
        lif.set("g_reset", neuron_config["g_reset"])
        lif.set("g_leak", neuron_config["g_leak"])
        lif.set("tau_refrac", neuron_config["tau_refrac"])
        decor.paint('(all)', lif)
        decor.place('"center"', arbor.spike_detector(v_thresh), "spike_detector")

        # plastic excitatory synapse
        syn_config_stdp = self.config["synapses"]["cond_exp_stdp"]
        mech_expsyn_exc = arbor.mechanism('expsyn_stdp')
        mech_expsyn_exc.set('tau', syn_config_stdp["tau"])
        mech_expsyn_exc.set('e', syn_config_stdp["reversal_potential"])
        mech_expsyn_exc.set('taupre', syn_config_stdp["tau_pre"])
        mech_expsyn_exc.set('taupost', syn_config_stdp["tau_post"])
        mech_expsyn_exc.set('Apre', syn_config_stdp["A_pre"])
        mech_expsyn_exc.set('Apost', syn_config_stdp["A_post"])
        mech_expsyn_exc.set('max_weight', 50)
        decor.place('"center"', mech_expsyn_exc, "expsyn_stdp_exc")

        # inhibitory synapse
        syn_config = self.config["synapses"]["cond_exp"]
        mech_expsyn_inh = arbor.mechanism('expsyn')
        mech_expsyn_inh.set('tau', syn_config["tau"])
        mech_expsyn_inh.set('e', syn_config["reversal_potential"])
        decor.place('"center"', mech_expsyn_inh, "expsyn_inh")

        return arbor.cable_cell(tree, labels, decor)

    def event_generators(self, gid):
        """Return event generators on gid."""
        assert gid == 0

        syn_config_stdp = self.config["synapses"]["cond_exp_stdp"]
        syn_config = self.config["synapses"]["cond_exp"]

        stimulus_times_exc = syn_config_stdp["stimulus_times"]
        stimulus_times_inh = syn_config["stimulus_times"]

        spike_exc = arbor.event_generator(
            "expsyn_stdp_exc",
            syn_config_stdp["weight"],
            arbor.explicit_schedule(stimulus_times_exc))
        spike_inh = arbor.event_generator(
            "expsyn_inh", syn_config["weight"], arbor.explicit_schedule(stimulus_times_inh))

        return [spike_exc, spike_inh]

    def probes(self, gid):
        """Return probes on gid."""
        assert gid == 0

        return [arbor.cable_probe_membrane_voltage('"center"'),
                arbor.cable_probe_point_state(0, "expsyn_stdp", "g"),
                arbor.cable_probe_point_state(1, "expsyn", "g"),
                arbor.cable_probe_point_state(0, "expsyn_stdp", "weight_plastic")]

    def global_properties(self, kind):
        """Return the global properties."""
        assert kind == arbor.cell_kind.cable

        return self.the_props


def main(variant):
    """Runs simulation and stores results."""

    # set up simulation and run
    config = json.load(open(f"config_{variant}.json", 'r'))
    recipe = SingleRecipe(config)

    context = arbor.context()
    domains = arbor.partition_load_balance(recipe, context)
    sim = arbor.simulation(recipe, domains, context)

    sim.record(arbor.spike_recording.all)

    reg_sched = arbor.regular_schedule(config["simulation"]["dt"])
    handle_mem = sim.sample((0, 0), reg_sched)
    handle_ge = sim.sample((0, 1), reg_sched)
    handle_gi = sim.sample((0, 2), reg_sched)
    handle_weight_plastic = sim.sample((0, 3), reg_sched)

    sim.run(tfinal=config["simulation"]["runtime"],
            dt=config["simulation"]["dt"])

    # readout traces and spikes
    data_mem, _ = sim.samples(handle_mem)[0]
    data_ge, _ = sim.samples(handle_ge)[0]
    data_gi, _ = sim.samples(handle_gi)[0]
    data_weight_plastic, _ = sim.samples(handle_weight_plastic)[0]

    # collect data and store
    data_stacked = numpy.column_stack(
        [data_mem[:, 0], data_mem[:, 1], data_ge[:, 1], data_gi[:, 1], data_weight_plastic[:, 1]])

    spike_times = sorted([s[1] for s in sim.spikes()])

    numpy.savetxt(f'arbor_traces_{variant}.dat', data_stacked)
    numpy.savetxt(f'arbor_spikes_{variant}.dat', spike_times)


if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('variant', help="name of variant, e.g., brian2_arbor")
    args = parser.parse_args()

    main(args.variant)
