Test data for the network synthesis process

../bin/network_synthesis

The test data consists of an impedance which may be in rational function form or alternatively
in pole-residue form (as produced by the Vector_fit process)

The script run_test runs the network synthesis process. The argument of run_test is the filename for 
the test data. The script runs the network synthesis process and produces a circuit file (ngspice_circuit.cir)
suitable for running in ngspice. The ngspice circuit calculates the frequeny domain voltage due to 
a unit current source and hence gives the (complex) impedance. 
A process 'compare_results' provides a measure of the difference between the complex spice model impedance and the
impedance calculated by evaluating the input rational function or pole-zero impedance function. 
Gnuplot then plots the comparison between the evaluation of the input function and the ngspice result. 

Examples:

1. Rational function describing the impedance of a ladder network with R in series with a parallel combination
of a C branch and a RL branch.

run_tet test_RCRL_ladder

