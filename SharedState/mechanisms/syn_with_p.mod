NEURON {
    POINT_PROCESS syn_with_p
    NONSPECIFIC_CURRENT i
	USEION p WRITE pi WRITE po
}

ASSIGNED {}

INITIAL {
	pi = 0
}

STATE {
    pi
}

BREAKPOINT {
    SOLVE state METHOD cnexp
}

DERIVATIVE state {
    pi' = 0
}

NET_RECEIVE(weight) {
    pi = pi + weight
}
