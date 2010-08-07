#!/usr/bin/gnuplot

# plot file for speedup
set data style lp
#set logscale x 2
set key top left
#set grid


#set ytics ("1x" 1, "2x" 2, "3x" 3, "4x" 4, "5x" 5, \
        #"6x" 6, "7x" 7, "8x" 8) 

set xlabel 'Number of cores' font "Helvetica,16"
set ylabel "GFLOPs" font "Helvetica,16"

set title "DGEMM perf - PLuTo vs. Rest of the Universe - Intel Core 2 Quad Q6600" font "Helvetica,20"

set xtics (1, 2, 3, 4)

set xrange [0.7:4.3]
set yrange [0:35]

set terminal postscript enhanced color eps #"Times-Roman" 22

set output 'matmul-par.eps'

set style line 1 lt 9 lw 3 pt 3 ps 0.5
set style line 2 lt 7 lw 2 pt 3 ps 0.5

plot 'matmul-par.dat' using 1:(2048/1000*2048/1000*2048/1000*2/$2) title 'PLuTo 0.0.1' ,  \
         '' using 1:(2048/1000*2048/1000*2048/1000*2/$3) title 'ICC -fast', \
         '' using 1:(2048/1000*2048/1000*2048/1000*2/$4) title 'ATLAS 3.8.0 (kSelMM)', \
         '' using 1:(2048/1000*2048/1000*2048/1000*2/$4*0.5222) title 'ATLAS 3.8.0 kGenMM', \
         '' using 1:(2048/1000*2048/1000*2048/1000*2/$5) title 'Intel MKL 9.1', \
         '' using 1:(2048/1000*2048/1000*2048/1000*2/$6) title 'NetLib BLAS 3.1 FC7' ls 2, \
         '' using 1:(2048/1000*2048/1000*2048/1000*2/$7) title 'Goto BLAS 1.23', \
         '' using 1:(2048/1000*2048/1000*2048/1000*2/$8) title 'PLuTo+ancc'
         #'' using 1:(4*9.6) title 'Absolute machine peak' ls 1
