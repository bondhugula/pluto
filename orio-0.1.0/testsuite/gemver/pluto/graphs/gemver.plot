#!/usr/bin/gnuplot

# plot file for speedup
set data style lp
#set logscale x 2
#set key top left
#set grid

#set ytics ("1x" 1, "2x" 2, "3x" 3, "4x" 4, "5x" 5, \
        #"6x" 6, "7x" 7, "8x" 8) 

set xlabel 'Matrix size (N)' #font "Helvetica,16"
set ylabel "GFLOPs" #font "Helvetica,16"

set title "GEMVER (Sequential)" #font "Helvetica,18"

set xtics ("2k" 2000, "4k" 4000, "6k" 6000, "8k" 8000, "10k" 10000)

set xrange [1000:11000]
set yrange [0:2]

set terminal postscript enhanced color eps #"Times-Roman" 22
#set terminal png enhanced tiny size 450,350

set output 'gemver.eps'
#set output 'gemver.png'

set style line 1 lt 9 lw 3 pt 3 ps 0.5
set style line 2 lt 7 lw 2 pt 3 ps 0.5

plot 'gemver.dat' using 1:2 title 'ICC -fast' ,  \
         '' using 1:3 title 'PLuTo 0.0.1', \
         '' using 1:4 title 'PLuTo+ancc'


