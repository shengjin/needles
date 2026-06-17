set term png
set out "a.png"
#set log y
#set log x
#set xrange [0:990000]
set yrange [400:500]
p "results/chain7.dat" u ($14) w p
