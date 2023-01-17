reset
set term gif animate delay 20
set size square

range = 100
set pm3d map
set cbrange [0.2:1.0]
set xrange [-range:range]
set yrange [-range:range]
set zrange [-range:range]
set xlabel "kpc"
set ylabel "kpc"
set ticslevel 0
set palette rgbformulae 31,13,10

file_name = "density_src.data"
last_chunk_index = 70
skip = 2
# 刻み幅の合計。xとyで同じ想定がされている。
step_num = 200

do for [i=0:last_chunk_index:skip]{
    start = i*(step_num*(step_num+1))+i*3+1

    time = system(sprintf('awk "NR == %d" %s',start,file_name))

    set title sprintf('t = %s',time)

    splot file_name using 1:2:(log10($3)) every :::1 index i with pm3d title "density"
}

# 入力ファイルのformat。異なる系列のデータの間は空白2行、同系列は1行

# t_n
#
# xn yn 値
# xn yn+1 値
#
# xn+1 yn 値
# xn+1 yn+1 値
# 
#
#t_n+1
#
# xn yn 値
# xn yn+1 値
#
# xn+1 yn 値
# xn+1 yn+1 値