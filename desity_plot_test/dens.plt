reset
set term gif animate
set size square
set output "sample.gif"

set pm3d map
set cbrange [-1:1]
set xrange [-1.5:1.5]
set yrange [-1.5:1.5]
set zrange [-1.5:1.5]
set ticslevel 0
set palette defined (-1 "blue", 0 "white", 1 "red")
splot "sample.data" with pm3d title "density"