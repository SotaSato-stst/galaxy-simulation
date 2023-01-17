reset
set term gif animate delay 10
set size square

set datafile separator ","
range = 100

set xrange [-range:range]
set yrange [-range:range]
set ticslevel 0

last_chunk_index = 100
skip = 1

do for [i=0:last_chunk_index:skip]{
    start = i*(num+3)+1
    stop  = (num*3)*(i+1) - 2

    time = system(sprintf('awk "NR == %d" %s',start,file_name))

    set title sprintf('t = %s',time)

    plot file_name every ::1 index i using 1:2 title "point" with points pt 7 lc 6 pointsize 0.2, factor_file_name index i using 1:2:3 title "point" with points pt 7 lc 2 pointsize 1.0
}

if(last_chunk_index%skip != 0){
    i     = last_chunk_index
    start = i*(num+3)+1
    stop  = (num*3)*(i+1) - 2

    time = system(sprintf('awk "NR == %d" %s',start,file_name))

    set title sprintf('t = %s',time)

    splot file_name every ::1 index i using 1:2:3 title "point" with points pt 7 lc 6 pointsize 0.2
}

print("")
