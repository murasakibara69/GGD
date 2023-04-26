set xlabel "x, м"
set ylabel "T, °C"
set nokey
set grid
set terminal pngcairo size 600, 400
set output "result.png"
plot "result.dat" u 1:2 w l
