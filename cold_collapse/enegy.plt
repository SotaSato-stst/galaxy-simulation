reset
set terminal png
set output "enegy.png"

set datafile separator ","
set xrange [0:5.0]
set yrange [-10000000000000:10000000000000]
set xlabel "t [Myr]"
set ylabel "Enegy [M_{solar} kpc^2 Myr^{-2}]"

file_name = "../cold_collapse/cold-collapse_enegy.csv"

plot file_name using 1:2 title " kinetic energy" with line,  file_name using 1:3 title "potential energy" with line, file_name using 1:4 title "mechanical energy" with line
