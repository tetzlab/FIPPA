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
        self.the_props.catalogue = self.the_cat

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
        stimulus_config = self.config["stimulus"]

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
        decor.place('"center"', arbor.threshold_detector(v_thresh), "spike_detector")

        decor.place('"center"', arbor.iclamp(stimulus_config["amplitude"]/1e-9), "iclamp")

        return arbor.cable_cell(tree, decor, labels)

    def event_generators(self, gid):
        """Return event generators on gid."""
        assert gid == 0

        return []

    def probes(self, gid):
        """Return probes on gid."""
        assert gid == 0

        return [arbor.cable_probe_membrane_voltage('"center"'),
                arbor.cable_probe_density_state('"center"', "adex", "w")
        ]

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
    sim = arbor.simulation(recipe, context, domains)

    sim.record(arbor.spike_recording.all)

    reg_sched = arbor.regular_schedule(config["simulation"]["dt"])
    handle_mem = sim.sample((0, 0), reg_sched)
    handle_w = sim.sample((0, 1), reg_sched)

    sim.run(tfinal=config["simulation"]["runtime"],
            dt=config["simulation"]["dt"])

    # readout traces and spikes
    data_mem, _ = sim.samples(handle_mem)[0]
    data_w, _ = sim.samples(handle_w)[0]

    # collect data and store
    data_stacked = numpy.column_stack(
        [data_mem[:, 0],
         numpy.clip(data_mem[:, 1], a_min=None, a_max=config["neuron"]["v_thresh"]),
         data_w[:, 1],
        ]
    )

    spike_times = sorted([s[1] for s in sim.spikes()])

    numpy.savetxt(f'arbor_traces_{variant}.dat', data_stacked)
    numpy.savetxt(f'arbor_spikes_{variant}.dat', spike_times)


if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('variant', help="name of variant, e.g., brian2_arbor")
    args = parser.parse_args()

    main(args.variant)
