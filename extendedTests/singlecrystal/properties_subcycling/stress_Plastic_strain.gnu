#set term x11
set term png size 640,480
set output 'stress_Plastic_strain.png'

set title  'ParaDiS Stress/Strain'
set xlabel 'Plastic Strain'  
set ylabel 'Stress'  
unset key

set format x '%0.1e'
set format y '%0.1e'

plot 'stress_Plastic_strain' u 1:2 w lp lc rgb '#ad0000' lt 1 lw 1 pt 7 ps 0.7

