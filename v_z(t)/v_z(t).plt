reset
set size square
set terminal png
set output "sample.png"

set xrange [0:2]
set yrange [0:10]

file_name = "density_src.data"

plot file_name title "v_z(t)"
