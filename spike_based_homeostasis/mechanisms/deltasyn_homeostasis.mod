NEURON {
    POINT_PROCESS deltasyn_homeostasis
    RANGE w_init, w_max, dw_plus, dw_minus, delta_factor
    NONSPECIFIC_CURRENT I
}

UNITS {
	(ms) = (milliseconds)
	(nA) = (nanoamp)
}

PARAMETER {
    w_init = 0 (nA) : initial weight
    w_max = 2 (nA) : maximum weight
    dw_plus = 0.05 (nA) : weight change on pre spike
    dw_minus = 0.05 (nA) : weight change on post spike
    delta_factor = 1 : correction factor for timestep with delta function
}

STATE {
    psc (nA) : variable for postsynaptic current
    w_plastic (nA) : variable for plastic weight
}

INITIAL {
    psc = 0
    w_plastic = w_init
}

BREAKPOINT {
    SOLVE state METHOD cnexp

    : set current depending on whether there is a spike or not
    if (psc > 0) {
		I = -psc
		psc = 0
	}
	else {
		I = 0
	}

    : apply upper weight bound
    if (w_plastic > w_max) {
        w_plastic = w_max
    } else if (w_plastic < 0) {
        w_plastic = 0
    }
}

DERIVATIVE state {
}

NET_RECEIVE(weight) {
    : apply postsynaptic current
    psc = w_plastic * delta_factor

    : update weight
    w_plastic = w_plastic + dw_plus
}

POST_EVENT(time) {
    : update weight
    w_plastic = w_plastic - dw_minus
}
