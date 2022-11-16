reset
set term gif animate
set output "keita.gif"

set datafile separator ","

set xrange [-1:1]
set yrange [-1:1]
set zrange [-1:1]
set view equal xyz
set ticslevel 0

num = 1000
last_chunk_index = 100
skip = 1

do for [i=0:last_chunk_index:skip]{
    start = i*(num+3)+1
    stop  = (num*3)*(i+1) - 2

    time = system(sprintf('awk "NR == %d" %s',start,"../plummer/plummer_keita.csv"))

    set title sprintf('t = %s',time)

    splot "../cold_collapse/cold-collapse.csv" every ::1 index i using 1:2:3 title "point" with points pt 7 lc 6 pointsize 1
}

if(last_chunk_index%skip != 0){
    i     = last_chunk_index
    start = i*(num+3)+1
    stop  = (num*3)*(i+1) - 2

    time = system(sprintf('awk "NR == %d" %s',start,"answer.csv"))

    set title sprintf('t = %s',time)

    splot "answer.csv" every ::1 index i using 1:2:3 title "point" with points pt 7 lc 6 pointsize 0.05
}

print("")