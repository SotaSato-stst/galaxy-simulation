reset
set term gif animate delay 10
set output "plummer_rho.gif"
set size square

set datafile separator ","
range = 3.0

set xrange [0:100]
set yrange [0:10000000000]
set ticslevel 0
set xlabel "x [kpc]"
set ylabel "density [M_{solar} kpc^{-3}]"
set logscale

num = 50
last_chunk_index = 340
skip = 1
file_name = "../plummer/rho.csv"

do for [i=0:last_chunk_index:skip]{
    start = i*(53)+1

    time = system(sprintf('awk "NR == %d" %s',start,file_name))

    set title sprintf('t = %s [Myr]',time)

    plot file_name every ::1 index i using 1:2 title "numerical solution" with line, (3 * (10 ** 9) / (4 * 3.14)) *  ((1+ (x) ** 2) ** (-2.5)) title "analytical solution"
}

print("")
