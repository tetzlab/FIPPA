: Simple spiker neuron
: SPikes when refactory period.
: The spike detector threshold should match the threshold of the mechanism.

NEURON {
    SUFFIX simple_spiker
    NONSPECIFIC_CURRENT i
    RANGE t_spike, tau_refrac, e_reset, e_thresh
}

UNITS {
    (mV) = (millivolt)
    (S) = (siemens)
}

STATE {
    refractory_counter : refractory time left
    t : current time step
}

INITIAL {
    refractory_counter = tau_refrac + 1 : start not refractory
    t = 0
}

PARAMETER {
    g = 1000 (S/cm2) : conductance (large, for fast dynamics)

    e_reset = 0 : reset potential
    e_thresh = 1 : spike threshold

    tau_refrac = 0.1 (ms) : refractory period

    t_spike = 0 : requested spike time
}

BREAKPOINT {
    SOLVE state METHOD cnexp
    LOCAL e

    : threshold crossed -> start refractory counter
    if (v > e_thresh) {
       refractory_counter = 0
    }

    : choose between spike and reset potential
    if (refractory_counter <= tau_refrac) {
       e = e_reset : reset by pulling potential back
    } else if (t > t_spike) {
       e = 1.1 * e_thresh : trigger spike by choosing potential just above threshold
    }

    i = -g*e
}

DERIVATIVE state {
    refractory_counter' = 1
    t' = 1
}
