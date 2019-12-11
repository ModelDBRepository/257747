from neuron import h, gui
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx

dtype = np.float64

# one-compartment cell (soma)
soma        = h.Section(name='soma')
soma.diam   = 50         # micron
soma.L      = 63.66198   # micron, so that area = 10000 micron2
soma.nseg   = 1          # adimensional
soma.cm     = 1          # uF/cm2
soma.Ra     = 70         # ohm-cm

soma.nseg   = 1
soma.insert('na15')      # insert mechanism
soma.ena    = 65
h.celsius   = 24         # temperature in celsius
v_init      = -120       # holding potential   
h.dt        = 0.01       # ms - value of the fundamental integration time step, dt, used by fadvance().

# clamping parameters
dur         = 500        # clamp duration, ms
step        = 3          # voltage clamp increment
st_cl       = -120       # clamp start, mV
end_cl      = 1          # clamp end, mV
v_cl        = -120       # actual voltage clamp, mV

#number of elements of the vector containing the values from st_cl to end_cl with the fixed step
L=len(np.arange(st_cl, end_cl, step))

# vectors for data handling
t_vec       = h.Vector() # vector for time
v_vec       = h.Vector() # vector for voltage
v_vec_t     = h.Vector() # vector for voltage as function of time
i_vec       = h.Vector() # vector for current 
ipeak_vec   = h.Vector() # vector for peak current
inorm_vec   = h.Vector() # vector for normalized current

# saving data (comment the following 4 lines if you don't want to save the data)
f1 = open('2_f_inact_v_vec.dat', 'w')
f2 = open('2_f_inact_inorm_vec.dat', 'w')
f1.write("voltage=[\n")
f2.write("normalized_current=[\n")


# a two-electrodes voltage clamp
f3cl = h.VClamp(soma(0.5))
f3cl.dur[0] = 40	     # ms
f3cl.amp[0] = -120	     # mV
f3cl.dur[1] = dur        # ms
f3cl.amp[1] = v_cl       # mV
f3cl.dur[2] = 20         # ms
f3cl.amp[2] = -10        # mV


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


#figure definition
fig = plt.figure(figsize=(20,15))
fig.suptitle('2. Fast inactivation availability', fontsize=15, fontweight='bold')

ax1 = plt.subplot2grid((2, 4), (0, 0), colspan=2)
ax1.set_xlim(0,560)
ax1.set_ylim(-121,10)
ax1.set_xlabel('Time $(ms)$')
ax1.set_ylabel('Voltage $(mV)$')
ax1.set_title('Time/Voltage relation')

ax2 = plt.subplot2grid((2,4), (0, 2))
ax2.set_xlim(538,548)
ax2.set_ylim(-1.5,0.1)
ax2.set_xlabel('Time $(ms)$')
ax2.set_ylabel('Current density $(mA/cm^2)$')
ax2.set_title('Time/Current density relation - zoom in')

ax3 = plt.subplot2grid((2,4), (0, 3))
ax3.set_xlim(538,548)
ax3.set_ylim(-0.015,0.001)
ax3.set_xlabel('Time $(ms)$')
ax3.set_ylabel('Current density $(mA/cm^2)$')
ax3.set_title('Time/Current density relation - zoom in')

ax4 = plt.subplot2grid((2,4), (1, 0), colspan=2)
ax4.set_xlim(0,560)
ax4.set_ylim(-1.5,0.25)
ax4.set_xlabel('Time $(ms)$')
ax4.set_ylabel('Current density $(mA/cm^2)$')
ax4.set_title('Time/Current density relation')

ax5 = plt.subplot2grid((2,4), (1, 2), colspan=2)
ax5.set_xlim(-125,3)
ax5.set_ylim(-0.05,1.05)        
ax5.set_xlabel('Voltage $(mV)$')
ax5.set_ylabel('Normalized current')
ax5.set_title('Voltage/Normalized current relation')



fig.subplots_adjust(wspace=0.5)
fig.subplots_adjust(hspace=0.5)

# to plot in rainbow colors
values=range(L)
rbw = cm = plt.get_cmap('rainbow') 
cNorm  = colors.Normalize(vmin=0, vmax=values[-1])
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=rbw)

# clamping definition
def clamp(v_cl):

    f3cl.amp[1] = v_cl
    h.finitialize(v_init)  # calling the INITIAL block of the mechanism inserted in the section.

    # parameters initialization
    peak_curr = 0
    dens = 0
    t_peak = 0
    
    while (h.t<h.tstop): # runs a single trace, calculates peak current
        dens = f3cl.i/soma(0.5).area()*100.0-soma(0.5).i_cap # clamping current in mA/cm2, for each dt
        t_vec.append(h.t)       # code for store the current
        v_vec_t.append(soma.v)  # trace to be plotted
        i_vec.append(dens)      # trace to be plotted
        
        if ((h.t>=540)and(h.t<=542)):     # evaluate the peak (I know it is there)
            if(abs(dens)>abs(peak_curr)):
                peak_curr = dens        
                t_peak = h.t
                
        h.fadvance()

    # updates the vectors at the end of the run        
    v_vec.append(v_cl)             
    ipeak_vec.append(peak_curr)


### start program

def start():
    h.tstop = 40 + dur + 20 
    v_vec.resize(0)
    ipeak_vec.resize(0)


    k=0     # counter    
    for v_cl in np.arange(st_cl, end_cl, step): # iterates across voltages

        print('Voltage Clamp:    ', v_cl,'mV')

        # resizing the vectors
        t_vec.resize(0)
        i_vec.resize(0)
        v_vec_t.resize(0) 
       

        clamp(v_cl)
            
        # code for showing traces
        colorVal1 = scalarMap.to_rgba(v_cl-st_cl-k*(step-1)) 
        k=k+1

        ln1,=ax1.plot(t_vec, v_vec_t,color=colorVal1)
        ln2,=ax2.plot(t_vec, i_vec,color=colorVal1)
        ln3,=ax3.plot(t_vec, i_vec,color=colorVal1)
        ln4,=ax4.plot(t_vec, i_vec,color=colorVal1)

    ipeak_min = ipeak_vec.min()           # normalization of peak current with respect to the min since the values are negative

    for i in range(0, len(ipeak_vec), 1):
         colorVal2 = scalarMap.to_rgba(i)
         inorm_vec.append(ipeak_vec.x[i]/ipeak_min)
         ln5,=ax5.plot(v_vec.x[i], inorm_vec.x[i], 'o', c=colorVal2)

         #printing and saving data (comment the following line if you don't want to print the data)
         print('Voltage:   ', v_vec.x[i],'mV', ',   Normalized current:   ', inorm_vec.x[i])
         # comment the following 2 lines if you don't want to save the data)
         f1.write("%s ,\n" % v_vec.x[i]) 
         f2.write("%s ,\n" % inorm_vec.x[i]) 

    #saving the figure (comment the following line if you don't want to save the figure)   
    plt.savefig('2. Fast inactivation availability', format='pdf', dpi=300, orientation='portrait')    

    # comment the following 4 lines if you don't want to save the data
    f1.write("];")
    f2.write("];")
    f1.close()
    f2.close()

    plt.show()


start()




