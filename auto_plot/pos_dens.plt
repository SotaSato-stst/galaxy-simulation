set terminal png

set pm3d map
set xrange [0:20000]
set yrange [-100:100]
set cbrange [0:3.0]
set palette rgbformulae 31,13,10

file_name = "density_src.data"

splot file_name using 1:2:(log10($3)) with pm3d title ""

