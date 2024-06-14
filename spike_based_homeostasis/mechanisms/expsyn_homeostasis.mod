NEURON {
    POINT_PROCESS expsyn_homeostasis
    RANGE tau, e, dw_plus, dw_minus, w_max
    NONSPECIFIC_CURRENT i
}

UNITS {
    (mV) = (millivolt)
}

PARAMETER {
    tau = 2.0 (ms) : synaptic time constant
    e = 0 (mV) : reversal potential
    dw_plus = 0.05 (uS) : weight change on pre spike
    dw_minus = 0.05 (uS) : weight change on post spike
    w_max = 2 (uS) : maximum weight
}

STATE {
    g
    weight_plastic
}

INITIAL {
    g = 0
    weight_plastic = 0
}

BREAKPOINT {
    SOLVE state METHOD cnexp
    i = g*(v-e)

    : apply weight bounds
    if (weight_plastic > w_max) {
        weight_plastic = w_max
    } else if (weight_plastic < 0) {
        weight_plastic = 0
    }
}

DERIVATIVE state {
    g' = -g/tau
}

NET_RECEIVE(weight) {
    g = max(0, weight + weight_plastic)
    weight_plastic = weight_plastic + dw_plus
}

POST_EVENT(time) {
    weight_plastic = weight_plastic - dw_minus
}
