reset
set terminal png
set output "virial.png"

set datafile separator ","
set xrange [0:350]
set yrange [0:1]
set xlabel "t [Myr]"
set ylabel "virial ratio"
file_name = "../plummer/enegy.csv"

plot file_name using 1:5 title "ratio" with line
