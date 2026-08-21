#set term x11
set term png size 640,480
set output 'time_Plastic_strain.png'

set title  'ParaDiS Simulation Time / Plastic Strain'
set xlabel 'Simulation Time'  
set ylabel 'Plastic Strain'  
unset key

set format x '%0.0e'
set format y '%0.0e'

plot 'time_Plastic_strain' u 1:2 w lp lc rgb '#ad0000' lt 1 lw 1 pt 7 ps 0.7

