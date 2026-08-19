#set term x11
set term png size 640,480
set output 'density.png'

set title  'ParaDiS Dislocation Density'
set xlabel 'Strain'  
set ylabel 'Dislocation Density'  
unset key

set format x '%0.0e'
set format y '%0.0e'

plot 'density' u 2:3 w lp lc rgb '#ad0000' lt 1 lw 1 pt 7 ps 0.7

