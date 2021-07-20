# Spike-based homeostasis

[![Spike-based homeostasis](https://github.com/tetzlab/FIPPA/actions/workflows/homeostasis.yml/badge.svg)](https://github.com/tetzlab/FIPPA/actions/workflows/homeostasis.yml)

With input from FIPPA, Arbor was extended with the possibility for event-driven plasticity.

This allows for spike-based homeostasis.

Following the approach in O. Breitwieser's thesis "Towards a
Neuromorphic Implementation of Spike-Based Expectation Maximization",
two poisson stimuli are connected to a neuron. One with a varying rate
and the other with a fixed rate.  The synaptic weight from the varying
rate stimulus to the neuron is fixed. The synaptic weight from the
fixed rate stimulus to the neuron is plastic and tries to keep the
neuron at a firing rate that is determined by the parameters of the
plasticity rule.

![Demonstration of spike-based homeostasis](homeostasis.svg)