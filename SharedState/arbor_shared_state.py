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

    def __init__(self):
        """Initialize the recipe from config."""

        # The base C++ class constructor must be called first, to ensure that
        # all memory in the C++ class is initialized correctly.
        arbor.recipe.__init__(self)

        self.the_props = arbor.neuron_cable_properties()
        self.the_cat = arbor.load_catalogue("./custom-catalogue.so")
        self.the_cat.extend(arbor.default_catalogue(), "")
        self.the_props.catalogue = self.the_cat

        self.the_props.set_ion('p', valence=0, int_con=0, ext_con=0, rev_pot=0)

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

        # morphology
        tree = arbor.segment_tree()
        radius = 1

        tree.append(arbor.mnpos,
                    arbor.mpoint(-radius, 0, 0, radius),
                    arbor.mpoint(radius, 0, 0, radius),
                    tag=1)

        labels = arbor.label_dict({'center': '(location 0 0.5)'})

        # cell mechanism
        decor = arbor.decor()

        pas_with_p = arbor.mechanism("pas_with_p")

        decor.paint('(all)', arbor.density(pas_with_p))

        decor.place('"center"', arbor.synapse("syn_with_p"), "syn_0")

        decor.place('"center"', arbor.synapse("syn_with_p"), "syn_1")

        return arbor.cable_cell(tree, labels, decor)

    def event_generators(self, gid):
        """Return event generators on gid."""
        assert gid == 0

        spike_times_0 = [5, 10, 15, 35]
        spike_times_1 = [20, 25, 30, 35]

        return [arbor.event_generator("syn_0", 1., arbor.explicit_schedule(spike_times_0)),
                arbor.event_generator("syn_1", 1., arbor.explicit_schedule(spike_times_1))
        ]

    def probes(self, gid):
        """Return probes on gid."""
        assert gid == 0

        return [arbor.cable_probe_membrane_voltage('"center"'),
                arbor.cable_probe_ion_int_concentration('"center"', "p"),
        ]

    def global_properties(self, kind):
        """Return the global properties."""
        assert kind == arbor.cell_kind.cable

        return self.the_props


def main():
    """Runs simulation and stores results."""

    # set up simulation and run
    recipe = SingleRecipe()

    context = arbor.context()
    domains = arbor.partition_load_balance(recipe, context)
    sim = arbor.simulation(recipe, domains, context)

    sim.record(arbor.spike_recording.all)

    dt = 0.1
    reg_sched = arbor.regular_schedule(dt)
    handle_mem = sim.sample((0, 0), reg_sched)
    handle_pi = sim.sample((0, 1), reg_sched)

    sim.run(tfinal=100, dt=dt)

    # readout traces and spikes
    data_mem, _ = sim.samples(handle_mem)[0]
    data_pi, _ = sim.samples(handle_pi)[0]

    # collect data and store
    data_stacked = numpy.column_stack(
        [data_mem[:, 0],
         data_mem[:, 1],
         data_pi[:, 1],
        ]
    )

    spike_times = sorted([s[1] for s in sim.spikes()])

    numpy.savetxt(f'arbor_traces.dat', data_stacked)
    numpy.savetxt(f'arbor_spikes.dat', spike_times)


if __name__ == '__main__':

    main()
