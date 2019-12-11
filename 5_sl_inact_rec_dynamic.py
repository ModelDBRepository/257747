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
soma.diam   = 50        # micron
soma.L      = 63.66198  # micron, so that area = 10000 micron2
soma.nseg   = 1         # dimensionless
soma.cm     = 1         # uF/cm2
soma.Ra     = 70        # ohm-cm

soma.nseg   = 1
soma.insert('na15')     # insert mechanism
soma.ena    = 65
h.celsius   = 24        # temperature in celsus
v_init      = -120      # holding potential
h.dt        = 0.075     # ms - value of the fundamental integration time step, dt, used by fadvance().

# clamping parameters
min_inter    = 0.1      # pre-stimulus starting interval
max_inter    = 10000    # pre-stimulus endinging interval
num_pts    	 = 50       # number of points in logaritmic scale
cond_st_dur  = 1000     # conditioning stimulus duration
res_pot		 = -120     # resting potential
dur          = 0.1

# vector containing 'num_pts' values equispaced between log10(min_inter) and log10(max_inter)
vec_pts = np.logspace(np.log10(min_inter), np.log10(max_inter), num=num_pts)
L = len(vec_pts)

# vectors for data handling
rec_vec      = h.Vector()   
time_vec     = h.Vector()   
log_time_vec = h.Vector()   
t_vec        = h.Vector()
v_vec_t      = h.Vector()
i_vec_t      = h.Vector()


# voltage clamp with "five" levels
f3cl = h.VClamp_plus(soma(0.5))
f3cl.dur[0] = 5	  		     # ms
f3cl.amp[0] = -120    		 # mV
f3cl.dur[1] = cond_st_dur    # ms
f3cl.amp[1] = -20   		 # mV
f3cl.dur[2] = dur     		 # ms
f3cl.amp[2] = res_pot		 # mV
f3cl.dur[3] = 20     		 # ms
f3cl.amp[3] = -20   		 # mV
f3cl.dur[4] = 5     		 # ms
f3cl.amp[4] = -120   		 # mV

# finding the "initial state variables values"
from state_variables import finding_state_variables
initial_values = [x for x in finding_state_variables(v_init,h.celsius)]

print('initial values=  ', initial_values)


for seg in soma:
    seg.na15.iC1=initial_values[0]
    seg.na15.iC2=initial_values[1]
    seg.na15.iO1=initial_values[2]
    seg.na15.iI1=initial_values[3]
    seg.na15.iI2=initial_values[4]

# figures definition
fig, ax = plt.subplots(2, 2,figsize=(18,6))  
ln0, = ax[0,0].plot([], [], '*')
ln1, = ax[0,1].plot([], [], '-')
ln2, = ax[1,0].plot([], [], '-')
ln3, = ax[1,1].plot([], [], '-')
fig.suptitle('Recovery from slow inactivation', fontsize=15, fontweight='bold')

fig.subplots_adjust(wspace=0.5)
fig.subplots_adjust(hspace=0.5)

def init():
    ax[0,0].set_xlim(-150, 5 + cond_st_dur + max_inter + 20 + 100)
    ax[0,0].set_ylim(-121,0)
    ax[0,0].set_xlabel('Time $(ms)$')
    ax[0,0].set_ylabel('Voltage $(mV)$')
    ax[0,0].set_title('Time/Voltage relation')

    ax[0,1].set_xlim(-50, 5 + cond_st_dur + max_inter + 20 + 100)
    ax[0,1].set_ylim(-1.75,0.2)
    ax[0,1].set_xlabel('Time $(ms)$')
    ax[0,1].set_ylabel('Current density $(mA/cm^2)$')
    ax[0,1].set_title('Time/Current relation')

    ax[1,0].set_xlim(-150, 5 + cond_st_dur + max_inter + 20 + 5)
    ax[1,0].set_ylim(-0.1, 1.1)
    ax[1,0].set_xlabel('Time $(ms)$')
    ax[1,0].set_ylabel('Peak current')
    ax[1,0].set_title('Time/Fractional recovery')

    ax[1,1].set_xlim(-1.1,4.1)
    ax[1,1].set_ylim(-0.1, 1.1)
    ax[1,1].set_xlabel('Log(Time)')
    ax[1,1].set_ylabel('Peak current')
    ax[1,1].set_title('Log(Time)/Fractional recovery')    


    return ln0, ln1, ln2, ln3,

#to plot in rainbow colors
values=range(L)
rbw = cm = plt.get_cmap('rainbow') 
cNorm  = colors.Normalize(vmin=0, vmax=values[-1])
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=rbw)


def animate(frame):
    dur2=vec_pts[int(frame)]

    #resizing vectors
    t_vec.resize(0) 
    i_vec_t.resize(0) 
    v_vec_t.resize(0) 

    
    Clamp(dur2)

    colorVal1 = scalarMap.to_rgba(values[int(frame)])
    colorVal2 = scalarMap.to_rgba(values[0:int(frame)+1])

    ln0,=ax[0,0].plot(t_vec, v_vec_t, color=colorVal1)
    ln1,=ax[0,1].plot(t_vec, i_vec_t, color=colorVal1)


    ln2=ax[1,0].scatter(time_vec, rec_vec, c=colorVal2)
    ln3=ax[1,1].scatter(log_time_vec, rec_vec, c=colorVal2)

    return ln0, ln1, ln2, ln3,


# clamping definition
def Clamp(dur): 

    f3cl.dur[2] = dur
    h.tstop = 5 + 1000 + dur + 20 + 5 
    h.finitialize(v_init)
   
    #variables initialization 
    pre_i1 = 0
    pre_i2 = 0
    dens   = 0        
     
    # initialization peak current
    peak_curr1 = 0    
    peak_curr2 = 0


    while (h.t<h.tstop): # runs a single trace, calculates peak current
        
        dens = f3cl.i/soma(0.5).area()*100.0-soma(0.5).i_cap # clamping current in mA/cm2, for each dt
        t_vec.append(h.t)      
        v_vec_t.append(soma.v) 
        i_vec_t.append(dens)


        if ((h.t>5)and(h.t<15)):       # evaluate the first peak
            if(pre_i1<abs(dens)):
                peak_curr1=abs(dens)
            pre_i1=abs(dens)

        
        if ((h.t>(5+cond_st_dur+dur))and(h.t<(15+cond_st_dur+dur))):  # evaluate the second peak 


            if(pre_i2<abs(dens)):
                peak_curr2=abs(dens)
            pre_i2=abs(dens)

        h.fadvance()

    

    if len(time_vec) > L-1:   # resizing vectors when the protocol is completed (it is needed for looping the animation)
        rec_vec.resize(0) 
        time_vec.resize(0)
        log_time_vec.resize(0)
        

    time_vec.append(dur)
    log_time_vec.append(np.log10(dur))
    rec_vec.append(peak_curr2/peak_curr1)

### start program

def start():
    k=0 #cunter

    for dur in vec_pts: 

        t_vec.resize(0)
        v_vec_t.resize(0) 
        rec_vec.resize(0) 
        time_vec.resize(0)
        log_time_vec.resize(0)
        k+=1

        # animtion     
        ani = animation.FuncAnimation(fig, animate, frames=L,
                         init_func=init, blit=True, interval=500, repeat=True)

    plt.show()


start()


















