from neuron import h, gui
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx

dtype = np.float64

# one-compartment cell (soma)
soma        = h.Section(name='soma')
soma.diam   = 50        # micron
soma.L      = 63.66198  # micron, so that area = 10000 micron2
soma.nseg   = 1         # adimensional
soma.cm     = 1         # uF/cm2
soma.Ra     = 70        # ohm-cm

soma.nseg   = 1
soma.insert('na15')     #insert mechanism
soma.ena    = 65
h.celsius   = 24        # temperature in celsius
v_init      = -120      # holding potential
h.dt        = 0.01      # ms - value of the fundamental integration time step, dt, used by fadvance().

# clamping parameters
dur         = 20        # clamp duration, ms
step        = 2         # voltage clamp increment, the user can 
st_cl       = -90       # clamp start, mV
end_cl      = 11        # clamp end, mV
v_cl        = -90       # actual voltage clamp, mV

#number of elements of the vector containing the values from st_cl to end_cl with the fixed step
L=len(np.arange(st_cl, end_cl, step))

# vectors for data handling
t_vec       = h.Vector()  # vector for time
v_vec       = h.Vector()  # vector for voltage
v_vec_t     = h.Vector()  # vector for voltage as function of time
i_vec       = h.Vector()  # vector for current 
ipeak_vec   = h.Vector()  # vector for peak current
gpeak_vec   = h.Vector()  # vector for peak conductance
gnorm_vec   = h.Vector()  # vector for normalized conductance

# saving data (comment the following 4 lines if you don't want to save the data)
f1          = open('1_V_G_relation_v_vec.dat', 'w')
f2          = open('1_V_G_relation_norm_conductance.dat', 'w')
f1.write("voltage=[\n")
f2.write("normalized conductance=[\n")


# a two-electrodes voltage clamp
f3cl = h.VClamp(soma(0.5))
f3cl.dur[0] = 5	      # ms
f3cl.amp[0] = -120	  # mV
f3cl.dur[1] = dur     # ms
f3cl.amp[1] = v_cl    # mV
f3cl.dur[2] = 5       # ms
f3cl.amp[2] = -120    # mV

# finding the "initial state variables values"
from state_variables import finding_state_variables
initial_values = [x for x in finding_state_variables(v_init,h.celsius)]

print('Initial values [C1, C2, O1, I1, I2]=  ', initial_values)


for seg in soma:
    seg.na15.iC1=initial_values[0]
    seg.na15.iC2=initial_values[1]
    seg.na15.iO1=initial_values[2]
    seg.na15.iI1=initial_values[3]
    seg.na15.iI2=initial_values[4]

# figure definition
fig, ax     = plt.subplots(1, 3,figsize=(18,6))  
fig.suptitle('1. Voltage-Normalized conductance relation', fontsize=15, fontweight='bold')
fig.subplots_adjust(wspace=0.5)

ax[0].set_xlim(0,30)
ax[0].set_ylim(-121,20)
ax[0].set_xlabel('Time $(ms)$')
ax[0].set_ylabel('Voltage $(mV)$')
ax[0].set_title('Time/Voltage relation')

ax[1].set_xlim(0,30)
ax[1].set_ylim(-1.5,0.5)
ax[1].set_xlabel('Time $(ms)$')
ax[1].set_ylabel('Current density $(mA/cm^2)$')
ax[1].set_title('Time/Current density relation')

ax[2].set_xlim(-100,20)
ax[2].set_ylim(-0.05,1.05)        
ax[2].set_xlabel('Voltage $(mV)$')
ax[2].set_ylabel('Normalized conductance')
ax[2].set_title('Voltage/Normalized conductance')

# to plot in rainbow colors
values      = range(L)
rbw         = cm = plt.get_cmap('rainbow') 
cNorm       = colors.Normalize(vmin=0, vmax=values[-1])
scalarMap   = cmx.ScalarMappable(norm=cNorm, cmap=rbw)

# clamping definition
def clamp(v_cl):
    
    curr_tr = 0  # initialization of peak current
    cond_tr = 0  # initialization of peak conductance

    h.finitialize(v_init) # calling the INITIAL block of the mechanism inserted in the section.

    # initialization of variables used to commute the peak current and conductance    
    pre_i  = 0     
    dens   = 0

    f3cl.amp[1] = v_cl    # mV
    


    while (h.t<h.tstop): # runs a single trace, calculates peak current
        dens = f3cl.i/soma(0.5).area()*100.0-soma(0.5).i_cap # clamping current in mA/cm2, for each dt
        
        t_vec.append(h.t)       # code for storing the current
        v_vec_t.append(soma.v)  # trace to be plotted
        i_vec.append(dens)      # trace to be plotted

        if ((h.t>5)and(h.t<=10)):       # evaluate the peak 
            if(abs(dens)>abs(pre_i)):
                cond_tr = soma.g_na15   # updates the peak conductance
                curr_tr = dens          # updates the peak current

                
        h.fadvance()
        pre_i = dens

    # updates the vectors at the end of the run        
    v_vec.append(v_cl)              
    gpeak_vec.append(cond_tr)
    ipeak_vec.append(curr_tr)


### start program

def start():

    h.tstop = 5 + dur + 5       # time stop

    # resizing the vectors
    v_vec.resize(0)          
    gpeak_vec.resize(0)
    ipeak_vec.resize(0)

    k=0     # counter
    for v_cl in np.arange(st_cl, end_cl, step): # iterates across voltages

            # resizing the vectors
            t_vec.resize(0)  
            i_vec.resize(0) 
            v_vec_t.resize(0) 

            print('Voltage Clamp:    ', v_cl,'mV')

            clamp(v_cl)
            colorVal1 = scalarMap.to_rgba(v_cl-st_cl-k*(step-1)) # rainbow printing setting
            k=k+1

            ln0,=ax[0].plot(t_vec, v_vec_t,color=colorVal1)
            ln1,=ax[1].plot(t_vec, i_vec,color=colorVal1)
            
 
    gpeak_max = gpeak_vec.max()           # maximum value of the conductance used to normalize the conductance vector
    for i in range(0, len(gpeak_vec), 1):

         colorVal2 = scalarMap.to_rgba(i) # rainbow printing setting

         gnorm_vec.append(gpeak_vec.x[i]/gpeak_max) # normalization of peak conductance
         ln2,=ax[2].plot(v_vec.x[i], gnorm_vec.x[i], 'o', c=colorVal2)

         # printing and saving data (comment the following line if you don't want to print the data)
         print('Voltage:   ', v_vec.x[i],'mV', ',   Normalized conductance:   ', gnorm_vec.x[i] )
         # comment the following 2 lines if you don't want to save the data)
         f1.write("%s ,\n" % v_vec.x[i])        
         f2.write("%s ,\n" % gnorm_vec.x[i])            

    #saving the figure (comment the following line if you don't want to save the figure)   
    plt.savefig('1. Voltage-Normalized conductance relation', format='pdf', dpi=300, orientation='portrait')    


    # comment the following 4 lines if you don't want to save the data
    f1.write("];")
    f2.write("];")
    f1.close()
    f2.close()

    plt.show()


start()

