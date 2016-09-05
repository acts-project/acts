set nokey
set size ratio -1

set title "Hits in barrel (|z| < 70 cm) in x-y plane"
set xlabel "x [cm]"
set ylabel "y [cm]"
set xrange [-70:70]
set yrange [-70:70]
set terminal pngcairo size 500,500 enhanced font 'Verdana,10'

set output "hits_xy.png"
plot "hits.txt" using 1:(abs($3) < 70 ? $2 : $2/0) with points pt 7 ps 0.1

set title "Hits in in r-z plane"
set xlabel "z [cm]"
set ylabel "r [cm]"
set xrange [-100:100]
set yrange [-70:70]
set terminal pngcairo size 714,500 enhanced font 'Verdana,10'
set output "hits_rz.png"
plot "hits.txt" using 3:(sgn($1) * sqrt($1*$1 + $2*$2)) with points pt 7 ps 0.1

set title "Hits in in Eta-Phi plane"
set size ratio -1
set xlabel "Eta"
set ylabel "Phi"
set xrange [-5:5]
set yrange [-pi:pi]
set terminal pngcairo size 714,500 enhanced font 'Verdana,10'
set output "hits_etaphi.png"
plot "hits.txt" using (-log(tan(atan2(sqrt($1*$1 + $2*$2),$3)/2))):(atan2($2,$1)) with points pt 7 ps 0.1
