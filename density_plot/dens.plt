reset
set term gif animate delay 100
set size square
set output "sample.gif"

set pm3d map
set cbrange [-1:1]
set xrange [-1.5:1.5]
set yrange [-1.5:1.5]
set zrange [-1.5:1.5]
set ticslevel 0
set palette defined (-1 "blue", 0 "white", 1 "red")

file_name = "sample.data"
last_chunk_index = 2
skip = 1
# 刻み幅の合計。xとyで同じ想定がされている。
step_num = 3

do for [i=0:last_chunk_index:skip]{
    start = i*(step_num*(step_num+1))+i*3+1

    time = system(sprintf('awk "NR == %d" %s',start,file_name))

    set title sprintf('t = %s',time)

    splot file_name every :::1 index i with pm3d title "density"
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