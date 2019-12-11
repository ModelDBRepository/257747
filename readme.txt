Readme file

The present files complement the paper:

- Andreozzi E, Carannante I, D'Addio G, Cesarelli M, Balbi P. Phenomenological models of NaV1.5.
A side by side, procedural, hands-on comparison between Hodgkin-Huxley and kinetic formalisms. Under submission

'Na15.mod' contains a five-states markovian kinetic model of Nav1.5 ionic voltage-gated channel.
'vclmp_pl.mod' contains a voltage clamp with "five" levels.
'state_variables.py' contains a function (finding_state_variables) which solves a linear system to find the initial values of the state variables.

The above-mentioned three files are necessary to run the following simulations (they have to be in the same folder as the stand-alone python files listed below) reproducing current-voltage relationship, fast inactivation availability, recovery from fast inactivation, development of slow inactivation and recovery from slow inactivation. 
  
1_I_V_relationship_static.py
1_I_V_relationship_dynamic.py
2_f_inact_avail_static.py
2_f_inact_avail_dynamic.py
3_f_inact_rec_static.py
3_f_inact_rec_dynamic.py
4_ons_sl_inact_static.py
4_ons_sl_inact_dynamic.py
5_sl_inact_rec_static.py
5_sl_inact_rec_dynamic.py

Each file is implemented in two modes: one static mode, and the corresponding dynamic mode which allows the user to better appreciate the electrophysiological protocols because they are shown in animation.

The 'static' files produce as output two .dat files containing the x and y values of the main plot e.g. voltage and normalized conductance for the 1_I_V_relationship_static.py. It also saves a .pdf file with a summary figure. Of course both options can be switched off by commenting the 
specified lines.

The python scripts are user-frendly, allowing the user to change the clamping parameters and other values such as holding potential (v_init) and temperature (h.celsius). Since the rate transitions depend on these last two values, their modification result in a modification in the state variables initial values. For this reason each file imports state_variables.py and calls the function finding_state_variables. 
To speed up the simulation it is possible to change h.dt which represents the integration time step (used by fadvance()). It is worth noticing that it is not possible to use variable step under h.VClamp mode. 





