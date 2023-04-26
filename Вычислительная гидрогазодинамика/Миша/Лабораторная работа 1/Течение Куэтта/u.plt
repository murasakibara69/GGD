set ylabel "y, m"
set nokey
set terminal pngcairo size 600, 400
set grid
set output "u.png"
plot "result.dat" u 1:2 w l
