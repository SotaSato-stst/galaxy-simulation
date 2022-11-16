reset
set term gif animate
set output "plummer_10000_v2.gif"

set datafile separator ","

set xrange [-100:100]
set yrange [-:100]
set zrange [-3:3]
set view equal xyz
set ticslevel 0

num = 300
last_chunk_index = 4000
skip = 50

do for [i=0:last_chunk_index:skip]{
    start = i*(num+3)+1
    stop  = (num*3)*(i+1) - 2

    time = system(sprintf('awk "NR == %d" %s',start,"../plummer/plummer_rho_10000.csv"))

    set title sprintf('t = %s',time)

    splot "../plummer/plummer_rho_10000.csv" every ::1 index i using 1:2:3 title "point" with points pt 7 lc 6 pointsize 0.2
}

if(last_chunk_index%skip != 0){
    i     = last_chunk_index
    start = i*(num+3)+1
    stop  = (num*3)*(i+1) - 2

    time = system(sprintf('awk "NR == %d" %s',start,"../plummer/plummer_animation_v3.csv"))

    set title sprintf('t = %s',time)

    splot "answer.csv" every ::1 index i using 1:2:3 title "point" with points pt 7 lc 6 pointsize 0.05
}

print("")