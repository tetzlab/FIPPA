: Delta synapse without plasticity

NEURON {
	POINT_PROCESS deltasyn
	RANGE psc_spike, delta_factor
	NONSPECIFIC_CURRENT I
}

UNITS {
	(ms) = (milliseconds)
	(nA) = (nanoamp)
}

PARAMETER {
	psc_spike = 1 (nA) : increase in postsynaptic current through one spike
    delta_factor = 1 : correction factor for timestep with delta function
}

ASSIGNED {
}

STATE {
	psc (nA) : instantaneous postsynaptic current (non-zero if there a spike in this timestep)
}

INITIAL {
	psc = 0
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	
	if (psc > 0) { : if there is a spike
		I = -psc
		psc = 0
	}
	else {
		I = 0
	}
}

DERIVATIVE state {
}

NET_RECEIVE(weight) {
	psc = psc_spike * delta_factor
}

