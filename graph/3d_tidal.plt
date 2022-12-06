reset
set term gif animate
set output "tidal.gif"

set datafile separator ","

set xrange [-300:300]
set yrange [-300:300]
set zrange [-300:300]
set view equal xyz
set ticslevel 0

num = 200
last_chunk_index = 50
skip = 1

do for [i=0:last_chunk_index:skip]{
    start = i*(num+3)+1
    stop  = (num*3)*(i+1) - 2

    time = system(sprintf('awk "NR == %d" %s',start,"../final/tidal.csv"))

    set title sprintf('t = %s',time)

    plot "../final/tidal.csv" every ::1 index i using 1:2 title "point" with points pt 7 lc 6 pointsize 0.05
}

if(last_chunk_index%skip != 0){
    i     = last_chunk_index
    start = i*(num+3)+1
    stop  = (num*3)*(i+1) - 2

    time = system(sprintf('awk "NR == %d" %s',start,"../final/tidal.csv"))

    set title sprintf('t = %s',time)

    splot "../final/tidal.csv" every ::1 index i using 1:2:3 title "point" with points pt 7 lc 6 pointsize 0.05
}

print("")
