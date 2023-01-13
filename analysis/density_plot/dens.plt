reset
set term gif animate delay 100
set size square
set output "sample.gif"

set pm3d map
set cbrange [0:1]
set xrange [-300:300]
set yrange [-300:300]
set zrange [-300:300]
set ticslevel 0
set palette defined (-1 '#ffffff', 0 '#000090',1 '#000fff',2 '#0090ff',3 '#0fffee',4 '#90ff70',5 '#ffee00',6 '#ff7000',7 '#ee0000',8 '#7f0000')

file_name = "density_src2.data"
last_chunk_index = 50
skip = 5
# 刻み幅の合計。xとyで同じ想定がされている。
step_num = 200

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