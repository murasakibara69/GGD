set xlabel "u, м/с"
set ylabel "y, м"
set nokey
set grid
set terminal pngcairo size 600, 400
set output "result.png"
plot "result.dat" u 1:2 w l 
