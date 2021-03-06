: Adaptive exponential leaky-integrate and fire mechanism
: When crossing the spike threshold, the membrane is pulled
: towards the reset potential with a (possibly) high conductance
: for the duration of the refactory period.
: The spike detector threshold should match the threshold of the mechanism.

NEURON {
    SUFFIX adex
    NONSPECIFIC_CURRENT i
    RANGE g_reset, g_leak, e_reset, e_leak, e_thresh, tau_refrac, delta_t, v_t, a, b, tau_w
}

UNITS {
    (mV) = (millivolt)
    (S) = (siemens)
}

STATE {
    refractory_counter
    w
    v_clipped
}

INITIAL {
    refractory_counter = tau_refrac + 1 : start not refractory
    w = 0
    v_clipped = v
}

PARAMETER {
    g_reset = 1000 (S/cm2) : conductance towards reset potential
    g_leak = 0.1 (S/cm2) : conductance towards leak potential

    e_reset = -58 (mV) : reset potential
    e_leak = -70  (mV) : leak potential
    e_thresh = 0 (mV) : spike threshold

    tau_refrac = 0.1 (ms) : refractory period

    delta_t = 2 (mV) : threshold slope factor
    v_t = -50 (mV) : effective threshold potential
    a = 0.016 (S/cm2) : adaptation voltage coupling
    b = 0.48 (mA/cm2) : spike-triggered adaptation
    tau_w = 300 (ms) : adaptation time constant
}

BREAKPOINT {
    SOLVE state METHOD cnexp

    LOCAL g
    LOCAL e
    LOCAL adex_current

    : threshold crossed -> start refractory counter
    if (v > e_thresh) {
       refractory_counter = 0
       w = w + b
    }

    : restrict membrane voltage to below threshold
    : otherwise w is spoiled
    v_clipped = min(e_thresh, v)

    : choose between leak and reset potential
    if (refractory_counter <= tau_refrac) {
       g = g_reset
       e = e_reset
       adex_current = 0
    } else {
       g = g_leak
       e = e_leak
       adex_current = - g_leak*delta_t*exp((v_clipped-v_t)/delta_t) + w
    }

    i = g*(v - e) + adex_current
}

DERIVATIVE state {
    refractory_counter' = 1
    w' = (a*(v_clipped-e_leak) - w)/tau_w
}
