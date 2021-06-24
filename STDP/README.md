# Event-driven plasticity

With input from FIPPA, Arbor was extended with the possibility for event-driven plasticity.

A synapse with spike-timing dependent plasticity (STDP) has been added to Arbor.

The implementation is validated against the Brian2 simulator.

# Howto

Adapt `run_arbor.sh` and `run_brian2.sh` to your environment.

```shell
make
```

Parameters can be changed in `config_brian2_arbor.json` or copied to a new configuration
file, e.g., `config_modified.json`.

Then call `make` like below:

```shell
make comparison_modified.png
```

![Comparison between Arbor and Brian2](comparison_brian2_arbor.png)

# Requirements

* [Brian2](https://briansimulator.org) >= 2.4.2
* [Arbor](https://github.com/arbor-sim/arbor) >= 0.5.2
* [LIF mechanism for Arbor](https://github.com/arbor-sim/arbor/pull/1517)
* [matplotlib](https://matplotlib.org) >= 3.4.1
