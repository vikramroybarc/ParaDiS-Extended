#set term x11
set term png size 640,480
set output 'stress_Plastic_strain.png'

set title  'ParaDiS Stress/Strain'
set xlabel 'Plastic Strain'  
set ylabel 'Stress'  
unset key


plot 'stress_Plastic_strain' u ($1*100):($2/1e6) w lp lc rgb '#ad0000' lt 1 lw 1 pt 7 ps 0.7 ,\
'stress_Plastic_strain_sma2' u ($1*100):($2/1e6) w lp lc rgb '#000000' lt 1 lw 2 pt 7 ps 0.7
