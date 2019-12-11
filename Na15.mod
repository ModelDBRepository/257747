TITLE Nav1.5 ionic voltage-gated channel with kinetic scheme

COMMENT
A five-state markovian kinetic model of ionic channel.
Part of a study on kinetic models.
Author: Piero Balbi, April 2019
ENDCOMMENT

NEURON {
	SUFFIX na15
	USEION na READ ena WRITE ina
	RANGE gbar, ina, g ,iC1,iC2,iO1,iI1,iI2
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}


PARAMETER {
	v (mV)
	ena (mV)
	celsius
	gbar  = 0.1	 (mho/cm2)

    iC1
    iC2
    iO1
    iI1
    iI2
	
	C1C2b2	  = 8
	C1C2v2    = -16
	C1C2k2	  = -9
	
	C2C1b1	  = 2
	C2C1v1    = -82
	C2C1k1	  = 5
	C2C1b2	  = 8
	C2C1v2    = -16
	C2C1k2	  = -9

	C2O1b2	  = 8
	C2O1v2    = -26
	C2O1k2	  = -9
	
	O1C2b1	  = 3
	O1C2v1    = -92
	O1C2k1	  = 5
	O1C2b2	  = 8
	O1C2v2    = -26
	O1C2k2	  = -9
	
	O1I1b1	  = 8
	O1I1v1	  = -50
	O1I1k1	  = 4
	O1I1b2	  = 6
	O1I1v2	  = 10
	O1I1k2	  = -100
	
	I1O1b1	  = 0.00001
	I1O1v1	  = -20
	I1O1k1	  = 10
	
	I1C1b1	  = 0.35
	I1C1v1	  = -122
	I1C1k1	  = 9
	
	C1I1b2	  = 0.04
	C1I1v2	  = -78
	C1I1k2	  = -10
	
	I1I2b2	  = 0.00018
	I1I2v2	  = -60
	I1I2k2	  = -5

	I2I1b1	  = 0.001825
	I2I1v1	  = -88
	I2I1k1	  = 31
	
}

ASSIGNED {
	ina  (mA/cm2)
	g   (mho/cm2)
	
	C1C2_a
	C2C1_a 
	C2O1_a  
	O1C2_a
	O1I1_a
	I1O1_a
	I1I2_a
	I2I1_a
	I1C1_a 
	C1I1_a
	
	Q10
}

STATE {
	C1
	C2
	O1
	I1
	I2
}


INITIAL {
	Q10 = 3^((celsius-20(degC))/10 (degC))
: valori delle variabili di stato calcolati da finding_state_variables
	C1 = iC1
	C2 = iC2
	O1 = iO1
	I1 = iI1
	I2 = iI2
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gbar * O1	: (mho/cm2)
	ina = g * (v - ena)   	: (mA/cm2)
}


FUNCTION rates2(v, b, vv, k) {
	rates2 = (b/(1+exp((v-vv)/k)))
}

DERIVATIVE states {
		rates(v)
        C1' = I1C1_a*I1+C2C1_a*C2-(C1I1_a+C1C2_a)*C1
        C2' = C1C2_a*C1+O1C2_a*O1-(C2C1_a+C2O1_a)*C2
        O1' = C2O1_a*C2+I1O1_a*I1-(O1C2_a+O1I1_a)*O1
        I1' = I2I1_a*(1 - (C1+C2+O1+I1)) + C1I1_a*C1+O1I1_a*O1-(I1C1_a+I1I2_a+I1O1_a)*I1
:       I2' = I1I2_a*I1-I2I1_a*I2
:		I2 = 1 - (C1+C2+O1+I1)
}


PROCEDURE rates(v(mV)) {
	C1C2_a = Q10*(rates2(v, C1C2b2, C1C2v2, C1C2k2))
	C2C1_a = Q10*(rates2(v, C2C1b1, C2C1v1, C2C1k1) + rates2(v, C2C1b2, C2C1v2, C2C1k2))
	C2O1_a = Q10*(rates2(v, C2O1b2, C2O1v2, C2O1k2))
	O1C2_a = Q10*(rates2(v, O1C2b1, O1C2v1, O1C2k1) + rates2(v, O1C2b2, O1C2v2, O1C2k2))
	O1I1_a = Q10*(rates2(v, O1I1b1, O1I1v1, O1I1k1) + rates2(v, O1I1b2, O1I1v2, O1I1k2))
	I1O1_a = Q10*(rates2(v, I1O1b1, I1O1v1, I1O1k1))
	I1C1_a = Q10*(rates2(v, I1C1b1, I1C1v1, I1C1k1))
	C1I1_a = Q10*(rates2(v, C1I1b2, C1I1v2, C1I1k2))
	I1I2_a = Q10*(rates2(v, I1I2b2, I1I2v2, I1I2k2))
	I2I1_a = Q10*(rates2(v, I2I1b1, I2I1v1, I2I1k1))
}



