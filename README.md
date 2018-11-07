# network_synthesis
Synthesis of spice models which reproduce specified s-domain impedance functions.

The inpout to the process is an s-domain impedance function which may be in rational function form (GGI_TLM_filter_fit 
from the GGI_TLM project (https://github.com/ggiemr/GGI_TLM)) or alternatively
in pole-residue form as produced by the Vector_fit process (https://github.com/chrissmartt/vector_fit)

The output is a Spice sub-circuit which reproduces the impedance function. 

If the input function does not represent a physical impedance  due to a negative resistance at some frequency
then a resistance is added in series such that the impedance function is physical (Re{Z}>0 for all frequencies))
and the process proceeds with this impedance. The Spice sub-circuit produced is for the corrected (physical) impedance.
However the Spice test circuit calculates the difference in voltage between the corrected impedance with a 1A source
and the added impedance with a 1A source so the voltage difference is that of the input function.

If the input function cannot be made physical by adding series resistance then an error occurs.
 
