set xlabel "Frequency (Hz)"
set ylabel "Function value"
set title "Compare rational function and continued fraction"
plot "Rational_funtion.fout" u 1:2 title "Re{H_rational}" w p,\
     "Rational_funtion.fout" u 1:3 title "Im{H_rational}" w p,\
     "Continued_fraction.fout" u 1:2 title "Re{H_CF}" w l,\
     "Continued_fraction.fout" u 1:3 title "Im{H_CF}" w l
pause -1

set title "Compare rational function and spice result fraction"
plot "Rational_funtion.fout" u 1:2 title "Re{H_rational}" w p,\
     "Rational_funtion.fout" u 1:3 title "Im{H_rational}" w p,\
     "spice_result.fout" u 2:3 title "Re{Spice_model_result}" w l,\
     "spice_result.fout" u 2:4 title "Im{Spice_model_result}" w l
pause -1
