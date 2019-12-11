# The following script contains a function 'finding_state_variables(v,celsius)', which takes as input # the holding potential in mV and the temperature in celsius and solves the ODE system, describing the # dynamic of the states, at equilibrium.  

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

dtype = np.float64


def finding_state_variables(v,celsius):

    global C1, C2, O1, I1, I2

    ### parameters
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

    ###
    def rates2(v,b,vv,k):
        return b/(1+np.exp((v-vv)/k))
    ###

    Q10 = 3**((celsius-20)/10)

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

    ### solving the ODE system at equilibrium taking into account the conservation law:
    ### C1+C2+O1+I1+I2 =1 ->
    ### solving the following linear system 

    A= np.array([[-(C1I1_a+C1C2_a),+C2C1_a,0,I1C1_a],
                 [C1C2_a,-(C2C1_a+C2O1_a),+O1C2_a,0],
                 [0,+ C2O1_a, -(O1C2_a+O1I1_a),+I1O1_a],
                 [(C1I1_a-I2I1_a), -I2I1_a, (O1I1_a-I2I1_a),-(I1C1_a+I1I2_a+I1O1_a+I2I1_a)]])


    b=np.array([0,0,0,-I2I1_a])


    x=np.linalg.solve(A,b)

    C1=x[0]
    C2=x[1]
    O1=x[2]
    I1=x[3]
    I2= 1-np.sum([x])

    return C1, C2, O1, I1, I2




