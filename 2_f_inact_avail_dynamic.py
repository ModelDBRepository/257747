from neuron import h, gui
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation  

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
t_vec       = h.Vector()     # vector for time
v_vec       = h.Vector()     # vector for voltage
v_vec_t     = h.Vector()     # vector for voltage as function of time
i_vec       = h.Vector()     # vector for current 
ipeak_vec   = h.Vector()     # vector for peak current

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

### definizione figure

fig, ax = plt.subplots(1,3,figsize=(18,4))  
ln0, = ax[0].plot([], [], '-')
ln1, = ax[1].plot([], [], '-')
ln2, = ax[2].plot([], [], '-')
fig.subplots_adjust(wspace=0.5)
fig.suptitle('2. Fast inactivation availability', fontsize=15, fontweight='bold')

def init():

    ax[0].set_xlim(0,560)
    ax[0].set_ylim(-121,10)
    ax[0].set_xlabel('Time $(ms)$')
    ax[0].set_ylabel('Voltage $(mV)$')
    ax[0].set_title('Time/Voltage relation')

    ax[1].set_xlim(0,560)
    ax[1].set_ylim(-1.5,0.25)
    ax[1].set_xlabel('Time $(ms)$')
    ax[1].set_ylabel('Current density $(mA/cm^2)$')
    ax[1].set_title('Time/Current density relation')

    ax[2].set_xlim(-121,10)
    ax[2].set_ylim(-1.5,0.25)
    ax[2].set_xlabel('Voltage $(mV)$')
    ax[2].set_ylabel('Current density $(mA/cm^2)$')
    ax[2].set_title('Current density/Voltage relation')



    return ln0, ln1, ln2, 

#to plot in rainbow colors
values    =range(L)
rbw       = cm = plt.get_cmap('rainbow') 
cNorm     = colors.Normalize(vmin=0, vmax=values[-1])
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=rbw)

# animation definition
def animate(frame):
        ran=np.arange(st_cl, end_cl, step)
        v_cl2=ran[int(frame)]

        #resizing vectors
        t_vec.resize(0) 
        i_vec.resize(0) 
        v_vec_t.resize(0) 

        clamp(v_cl2)

        colorVal1 = scalarMap.to_rgba(values[int(frame)])

        colorVal2 = scalarMap.to_rgba(values[0:int(frame)+2])
        
        ln0,=ax[0].plot(t_vec, v_vec_t,color=colorVal1)
        ln1,=ax[1].plot(t_vec, i_vec,color=colorVal1)


        ln2=ax[2].scatter(v_vec, ipeak_vec, c=colorVal2)
        return ln0, ln1, ln2


# clamping definition
def clamp(v_cl):

    f3cl.dur[1]=dur     # ms
    f3cl.amp[1]=v_cl    # mV
    
    h.finitialize(v_init)
    
    peak_curr = 0
    dens = 0
    t_peak = 0
    

    while (h.t<h.tstop): # runs a single trace, calculates peak current
        dens = f3cl.i/soma(0.5).area()*100.0-soma(0.5).i_cap # clamping current in mA/cm2, for each dt

        t_vec.append(h.t)               # code for store the current
        v_vec_t.append(soma.v)          # trace to be plotted
        i_vec.append(dens)              # trace to be plotted
        
        if ((h.t>=540)and(h.t<=542)):       # evaluate the peak
            if(abs(dens)>abs(peak_curr)):
                peak_curr = dens            # updates the peak current
                t_peak = h.t
                
        h.fadvance()



    if len(v_vec) > L-1:    #resizing v_vec and ipeak_vec when the protocol is completed (it is needed for looping the animation)
        v_vec.resize(0) 
        ipeak_vec.resize(0)

    v_vec.append(v_cl)              # updates the vectors at the end of the run
    ipeak_vec.append(peak_curr)

def start():
    
    h.tstop = 40 + dur + 20     # time stop

    v_vec.resize(0)
    ipeak_vec.resize(0)

    #animation
    ani = animation.FuncAnimation(fig, animate, frames=L,
                     init_func=init, blit=True, interval=500, repeat=True)

    
    plt.show()


start()






