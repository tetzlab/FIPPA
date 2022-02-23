NEURON {
    SUFFIX pas_with_p
    NONSPECIFIC_CURRENT i
    RANGE g
    GLOBAL e
	USEION p WRITE pi
}

UNITS {
    (mV) = (millivolt)
    (S) = (siemens)
}

STATE {
   pi
}

INITIAL {
   pi = 0
}

DERIVATIVE state {
   pi' = -pi/2.
}

PARAMETER {
    g = .001 (S/cm2)
    e = -70  (mV) : Taken from nrn
}

BREAKPOINT {
}

