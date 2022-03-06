#! /usr/bin/gnuplot

set term qt font "arial,12" enhanced #persist
set encoding utf8
set tmargin 1
set tics out
set grid
set nokey

plot "./prob.dat" using 1:2 with points pt 1 lt rgb "#FF0000"

pause mouse close

