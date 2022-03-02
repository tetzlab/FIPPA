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

    def __init__(self, config_soma, config_proximal, config_distal):
        """Initialize the recipe from config."""

        # The base C++ class constructor must be called first, to ensure that
        # all memory in the C++ class is initialized correctly.
        arbor.recipe.__init__(self)

        self.the_props = arbor.neuron_cable_properties()
        self.the_props.set_property(membrane_max=1)

        self.the_cat = arbor.load_catalogue("./custom-catalogue.so")
        self.the_cat.extend(arbor.default_catalogue(), "")
        self.the_props.catalogue = self.the_cat

        self.configs = {0 : config_soma,
                        1 : config_proximal,
                        2 : config_distal}

    def num_cells(self):
        """Return the number of cells."""
        return 3

    def num_sources(self, gid):
        """Return the number of spikes sources on gid."""
        assert gid in [0, 1, 2]
        return 1

    def num_targets(self, gid):
        """Return the number of post-synaptic targets on gid."""
        assert gid in [0, 1, 2]
        return 2

    def cell_kind(self, gid):
        """Return type of cell with gid."""
        assert gid in [0, 1, 2]
        return arbor.cell_kind.cable

    def cell_description(self, gid):
        """Return cell description of gid."""
        assert gid in [0, 1, 2]

        neuron_config = self.configs[gid]["neuron"]
        stimulus_config = self.configs[gid]["stimulus"]

        # morphology
        tree = arbor.segment_tree()
        radius = neuron_config["radius"]

        tree.append(arbor.mnpos,
                    arbor.mpoint(-radius, 0, 0, radius),
                    arbor.mpoint(radius, 0, 0, radius),
                    tag=1)

        surface = 4 * numpy.pi * radius**2

        labels = arbor.label_dict({'center': '(location 0 0.5)'})

        # cell mechanism
        decor = arbor.decor()

        cm = neuron_config["capacitance"] / (surface*1e-6**2)

        decor.set_property(Vm=neuron_config["e_leak"], cm=cm)
        adex = arbor.mechanism(neuron_config["type"])
        v_thresh = neuron_config["v_thresh"]

        # convert the leak/reset conductance from S to S/cm2 (surface is in um2)
        g_leak = neuron_config["g_leak"] / (surface*(1e-6**2)/(1e-2**2))

        g_reset = neuron_config["g_reset"] / (surface*(1e-6**2)/(1e-2**2))

        adex.set("e_thresh", v_thresh)
        adex.set("e_reset", neuron_config["v_reset"])

        adex.set("g_reset", g_reset)
        adex.set("g_leak", g_leak)

        adex.set("e_leak", neuron_config["e_leak"])

        adex.set("tau_refrac", neuron_config["tau_refrac"])
        adex.set("holdoff", neuron_config["holdoff"])

        adex.set("delta_t", neuron_config["delta_t"])
        adex.set("v_t", neuron_config["v_t"])

        # to S/cm2
        a = neuron_config["a"] / (surface*(1e-6**2)/(1e-2**2))
        adex.set("a", a)

        # to mA/cm2
        b = neuron_config["b"] / 1e-3 / (surface*(1e-6**2)/(1e-2**2))
        adex.set("b", b)

        adex.set("tau_w", neuron_config["tau_w"])

        decor.paint('(all)', arbor.density(adex))
        decor.place('"center"', arbor.spike_detector(v_thresh), "spike_detector")

        decor.place('"center"', arbor.iclamp(tstart=0, duration=stimulus_config["duration"], current=stimulus_config["amplitude"]/1e-9), "iclamp")

        junction_mech = arbor.junction('gj', {"g" : 0.1})
        decor.place('"center"', junction_mech, 'gj_label')

        return arbor.cable_cell(tree, labels, decor)

    def event_generators(self, gid):
        """Return event generators on gid."""
        assert gid in [0, 1, 2]

        return []

    def probes(self, gid):
        """Return probes on gid."""
        assert gid in [0, 1, 2]

        return [arbor.cable_probe_membrane_voltage('"center"')]

    def global_properties(self, kind):
        """Return the global properties."""
        assert kind == arbor.cell_kind.cable

        return self.the_props

    def gap_junctions_on(self, gid):
        assert gid in [0, 1, 2]

        connections = {0 : [1],
                       1 : [0, 2],
                       2 : [1]}

        return [arbor.gap_junction_connection((connection, 'gj_label'), 'gj_label', 1) for connection in connections[gid]]


def main(soma, proximal, distal):
    """Runs simulation and stores results."""

    # set up simulation and run
    configs = [json.load(open(f"config_{variant}.json", 'r')) for variant in [soma, proximal, distal]]
    recipe = SingleRecipe(*configs)

    context = arbor.context()
    domains = arbor.partition_load_balance(recipe, context)
    sim = arbor.simulation(recipe, domains, context)

    sim.record(arbor.spike_recording.all)

    assert len(set([c["simulation"]["runtime"] for c in configs])) == 1
    assert len(set([c["simulation"]["dt"] for c in configs])) == 1

    reg_sched = arbor.regular_schedule(configs[0]["simulation"]["dt"])
    handle_mem_soma = sim.sample((0, 0), reg_sched)
    handle_mem_proximal = sim.sample((1, 0), reg_sched)
    handle_mem_distal = sim.sample((2, 0), reg_sched)

    sim.run(tfinal=configs[0]["simulation"]["runtime"],
            dt=configs[0]["simulation"]["dt"])

    # readout traces and spikes
    data_mem_soma, _ = sim.samples(handle_mem_soma)[0]
    data_mem_proximal, _ = sim.samples(handle_mem_proximal)[0]
    data_mem_distal, _ = sim.samples(handle_mem_distal)[0]

    # collect data and store
    data_stacked = numpy.column_stack(
        [data_mem_soma[:, 0],
         data_mem_soma[:, 1],
         data_mem_proximal[:, 1],
         data_mem_distal[:, 1],
        ]
    )

    spike_times = sorted([s[1] for s in sim.spikes()])

    numpy.savetxt(f'arbor_traces_bac.dat', data_stacked)
    numpy.savetxt(f'arbor_spikes_bac.dat', spike_times)


if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('variants', help="soma, proximal and distal configuration", nargs=3)
    args = parser.parse_args()

    main(*args.variants)
