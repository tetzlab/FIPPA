#!/usr/bin/env python3.10
"""Arbor simulation of a single cell

Arbor simulation of a single cell receiving inhibitory and plastic
excitatory stimulus.
"""

import itertools
import json
import arbor
from arbor import units as U
import numpy

def generate_poisson_spike_train(rate, duration):
    """Returns a Poisson spike train.

    rate -- the rate of the spike train, in Hz
    duration -- the duration of the spike train, in s
    rng -- random number generator
    """

    n_events = numpy.round(rate * duration)

    if n_events:
        size = numpy.random.poisson(n_events)
        return numpy.sort(numpy.random.uniform(0.0, duration, size)) * U.s
    return numpy.array([]) * U.s

class SingleRecipe(arbor.recipe):
    """Implementation of Arbor simulation recipe."""

    def __init__(self, config, cat_file, considered_neuron):
        """Initialize the recipe from config."""

        # The base C++ class constructor must be called first, to ensure that
        # all memory in the C++ class is initialized correctly.
        arbor.recipe.__init__(self)

        # load custom catalogue
        catalogue = arbor.load_catalogue(cat_file)

        self.the_props = arbor.neuron_cable_properties()
        self.the_cat = catalogue
        self.the_cat.extend(arbor.default_catalogue(), "")
        self.the_props.catalogue = self.the_cat

        self.config = config
        self.considered_neuron = considered_neuron

    def num_cells(self):
        """Return the number of cells."""
        return 2

    def num_sources(self, gid):
        """Return the number of spikes sources on gid."""
        return 1

    def num_targets(self, gid):
        """Return the number of post-synaptic targets on gid."""
        return 2

    def cell_kind(self, gid):
        """Return type of cell with gid."""
        return arbor.cell_kind.cable

    def cell_description(self, gid):
        """Return cell description of gid."""

        sim_config = self.config["simulation"]
        neuron_config = self.config["neuron"]

        # morphology
        tree = arbor.segment_tree()
        radius = neuron_config["radius"]

        tree.append(arbor.mnpos,
                    arbor.mpoint(-radius, 0, 0, radius),
                    arbor.mpoint(radius, 0, 0, radius),
                    tag=1)

        labels = arbor.label_dict({'center': '(location 0 0.5)'})

        # current to current density conversion (necessary to simulate point neurons))
        height = 2*radius
        area = 2 * numpy.pi * radius*U.um * height*U.um # surface area of the cylinder (excluding the circle-shaped ends, since Arbor does not consider current flux there)
        area_cm2 = 2 * numpy.pi * radius * height * (1e-4)**2  # surface area of the cylinder in cm^2
        i_factor = (1e-9/1e-3) / area_cm2 # current to current density conversion factor (nA to mA/cm^2; necessary for point neurons)

        # cell mechanism
        dt = sim_config["dt"]*U.ms # timestep
        tau_mem = neuron_config["tau_e"]*U.ms # membrane time constant
        R_leak = 1/(1000*neuron_config["g_leak"])*U.MOhm # subthreshold membrane conductance
        C_mem = tau_mem/R_leak # absolute membrane capacitance for point neuron
        c_mem = C_mem / area # specific membrane capacitance
        v_thresh = neuron_config["v_thresh"]*U.mV # threshold voltage
        decor = arbor.decor()
        decor.set_property(Vm=neuron_config["e_leak"]*U.mV, cm=c_mem)
        lif = arbor.mechanism(neuron_config["type"])
        lif.set("e_thresh", v_thresh.value_as(U.mV))
        lif.set("e_leak", neuron_config["e_leak"])
        lif.set("e_reset", neuron_config["e_reset"])
        lif.set("g_leak", 1/R_leak.value_as(U.MOhm))
        lif.set("g_reset", neuron_config["g_reset"])
        lif.set("tau_refrac", neuron_config["t_ref"])
        lif.set("i_factor", i_factor)
        decor.paint('(all)', arbor.density(lif))
        decor.place('"center"', arbor.threshold_detector(v_thresh), "spike_detector")

        # plastic synapse with steady input
        syn_config_steady = self.config["stimulus"]["steady"]
        mech_steady = arbor.mechanism('deltasyn_homeostasis')
        mech_steady.set('dw_plus', syn_config_steady["dw_plus"])
        mech_steady.set('dw_minus', syn_config_steady["dw_minus"])
        mech_steady.set('w_init', syn_config_steady["w_init"])
        mech_steady.set('w_max', syn_config_steady["w_max"])
        mech_steady.set('delta_factor', tau_mem.value_as(U.ms)/dt.value_as(U.ms))
        decor.place('"center"', arbor.synapse(mech_steady), "input_steady")

        # static synapse with varying input
        syn_config_varying = self.config["stimulus"]["varying"]
        mech_varying = arbor.mechanism('deltasyn')
        mech_varying.set('psc_spike', syn_config_varying["w_non_plastic"])
        mech_varying.set('delta_factor', tau_mem.value_as(U.ms)/dt.value_as(U.ms))
        decor.place('"center"', arbor.synapse(mech_varying), "input_varying")

        return arbor.cable_cell(tree, decor, labels)

    def event_generators(self, gid):
        """Return event generators on gid."""

        syn_config_steady = self.config["stimulus"]["steady"]
        syn_config_varying = self.config["stimulus"]["varying"]
        dt_steady = syn_config_steady["dt"]*U.ms
        dt_varying = syn_config_varying["dt"]*U.ms

        # generate spikes for steady input
        stimulus_times_steady = numpy.concatenate([generate_poisson_spike_train(rate, dt_steady.value_as(U.s)) + dt_steady*i for i, rate in enumerate(syn_config_steady["rates"])])
        #numpy.savetxt(f"arbor_spikes_steady_{self.considered_neuron}.dat", stimulus_times_steady)

        # generate spikes for varying input
        stimulus_times_varying = numpy.concatenate([generate_poisson_spike_train(rate, dt_varying.value_as(U.s)) + dt_varying*i for i, rate in enumerate(syn_config_varying["rates"])])
        #numpy.savetxt(f"arbor_spikes_varying_{self.considered_neuron}.dat", stimulus_times_varying)
        
        # store varying input rate
        input_dts = numpy.arange(0., len(syn_config_varying["rates"])*dt_varying.value_as(U.ms), dt_varying.value_as(U.ms))
        input_x = list(itertools.chain(*zip(input_dts, input_dts)))
        input_y = [0] + list(itertools.chain(*zip(syn_config_varying["rates"], syn_config_varying["rates"])))[:-1]
        input_stacked = numpy.column_stack([input_x, input_y])
        numpy.savetxt(f"arbor_input_{self.considered_neuron}.dat", input_stacked)

        # create spike generators
        # neuron with homeostasis
        if self.considered_neuron == 1:
            spike_steady = arbor.event_generator("input_steady", 0, arbor.explicit_schedule(stimulus_times_steady))
        # neuron without homeostasis (control case)
        else:
            spike_steady = arbor.event_generator("input_steady", 0, arbor.explicit_schedule([]))
        spike_varying = arbor.event_generator("input_varying", 0, arbor.explicit_schedule(stimulus_times_varying)) # weight is set via 'psc_spike'

        return [spike_steady, spike_varying]

    def probes(self, gid):
        """Return probes on gid."""

        return [arbor.cable_probe_membrane_voltage('"center"', "tag_v"),
                arbor.cable_probe_point_state(0, "deltasyn_homeostasis", "psc", "tag_psc_0"),
                arbor.cable_probe_point_state(0, "deltasyn_homeostasis", "w_plastic", "tag_w_0"),
                arbor.cable_probe_point_state(1, "deltasyn", "psc", "tag_psc_1")]

    def global_properties(self, kind):
        """Return the global properties."""
        assert kind == arbor.cell_kind.cable

        return self.the_props


def main(config_file, catalogue, considered_neuron = 1):
    """Runs simulation and stores results.
    
    Parameters
    ----------
    config_file
        Name of the JSON configuration file.
    catalogue
        Name of the custom Arbor catalogue to use.
    considered_neuron
        The neuron to monitor: 0 for neuron without homeostasis, 1 for neuron with
        homeostasis.
    """

    # set up the simulation
    config = json.load(open(config_file, 'r'))
    recipe = SingleRecipe(config, catalogue, considered_neuron)

    context = arbor.context()
    domains = arbor.partition_load_balance(recipe, context)
    sim = arbor.simulation(recipe, context, domains)

    sim.record(arbor.spike_recording.all)

    reg_sched = arbor.regular_schedule(config["simulation"]["dt"]*U.ms)

    # probes on gid = 0 (neuron with homeostasis)
    handle_mem_with = sim.sample(0, "tag_v", reg_sched)
    handle_g_steady_with = sim.sample((0, "tag_psc_0"), reg_sched)
    handle_w_plastic_with = sim.sample((0, "tag_w_0"), reg_sched)
    handle_g_varying_with = sim.sample((0, "tag_psc_1"), reg_sched)

    # probes on gid = 1 (neuron without homeostasis)
    handle_mem_without = sim.sample((1, "tag_v"), reg_sched)
    handle_g_steady_without = sim.sample((1, "tag_psc_0"), reg_sched)
    handle_w_plastic_without = sim.sample((1, "tag_w_0"), reg_sched)
    handle_g_varying_without = sim.sample((1, "tag_psc_1"), reg_sched)

    # run the simulation
    sim.run(tfinal=config["simulation"]["runtime"]*U.ms,
            dt=config["simulation"]["dt"]*U.ms)

    # readout traces and spikes of neuron with homeostasis
    if considered_neuron == 1:
        data_mem, _ = sim.samples(handle_mem_with)[0]
        data_g_steady, _ = sim.samples(handle_g_steady_with)[0]
        data_w_plastic, _ = sim.samples(handle_w_plastic_with)[0]
        data_g_varying, _ = sim.samples(handle_g_varying_with)[0]
    # readout traces and spikes of neuron without homeostasis (control case)
    else:
        data_mem, _ = sim.samples(handle_mem_without)[0]
        data_g_steady, _ = sim.samples(handle_g_steady_without)[0]
        data_w_plastic, _ = sim.samples(handle_w_plastic_without)[0]
        data_g_varying, _ = sim.samples(handle_g_varying_without)[0]

    # collect data
    data_stacked = numpy.column_stack(
        [data_mem[:, 0], data_mem[:, 1], data_g_steady[:, 1], data_g_varying[:, 1], data_w_plastic[:, 1]])
    spike_times = []
    for s in sim.spikes():
        if s['source']['gid'] == considered_neuron:
            spike_times.append(s['time'])

    return data_stacked, spike_times


if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('config', help="name of config file")
    parser.add_argument('num_trials', type=int, help="number of trials to consider")
    parser.add_argument('considered_neuron', type=int, help="the type of neuron to consider")
    parser.add_argument('--catalogue', help="name of catalogue file library", default="homeostasis-catalogue.so")
    args = parser.parse_args()

    num_trials = args.num_trials
    data_stacked_sum = numpy.array([])
    spike_times_all = []
    for i in range(num_trials):
        data_stacked, spike_times = main(args.config, args.catalogue, args.considered_neuron)
        if data_stacked_sum.size == 0:
            data_stacked_sum = numpy.array(data_stacked)
        else:
            data_stacked_sum += data_stacked
        spike_times_all.extend(spike_times)
    
    numpy.savetxt(f'arbor_traces_{args.considered_neuron}.dat', data_stacked_sum / num_trials)
    numpy.savetxt(f'arbor_spikes_{args.considered_neuron}.dat', spike_times_all)
