reset
set term gif animate delay 10
set output "plummer.gif"
set size square

set datafile separator ","
range = 3.0

set xrange [-range:range]
set yrange [-range:range]
set ticslevel 0
set xlabel "x [kpc]"
set ylabel "y [kpc]"

num = 1000
last_chunk_index = 340
skip = 1
file_name = "../plummer/x_v.csv"

do for [i=0:last_chunk_index:skip]{
    start = i*(num+3)+1
    stop  = (num*3)*(i+1) - 2

    time = system(sprintf('awk "NR == %d" %s',start,file_name))

    set title sprintf('t = %s [Myr]',time)

    plot file_name every ::1 index i using 1:2 title "" with points pt 7 lc 6 pointsize 0.2
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
