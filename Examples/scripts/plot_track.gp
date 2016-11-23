set nokey
set size ratio -1

set title "Track momentum"
set xlabel "px [GeV]"
set ylabel "py [GeV]"
set xrange [-1.2:1.2]
set yrange [-1.2:1.2]
set terminal pngcairo size 500,500 enhanced font 'Verdana,10'
set output "track_pxpy.png"
plot "track.txt" using (cos($5)*$4):(sin($5)*$4)

set title "Track position"
set xlabel "x [mm]"
set ylabel "y [mm]"
set xrange [-1700:1700]
set yrange [-3400:0]
set terminal pngcairo size 714,500 enhanced font 'Verdana,10'
set output "track_xy.png"
plot "track.txt" using 1:2
