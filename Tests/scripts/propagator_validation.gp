set ylabel "Tracks"
set grid
set terminal pngcairo size 500,500 enhanced font 'Verdana,10'

unset logscale y
binwidth=0.001
bin(x,width)=width*floor(x/width)

set table 'tmp_x'
plot 'gsl_stepper_validation.txt' using (bin($1-$12,binwidth)):(1.0) smooth freq with boxes 

set table 'tmp_y'
plot 'gsl_stepper_validation.txt' using (bin($2-$13,binwidth)):(1.0) smooth freq with boxes 

set table 'tmp_z'
plot 'gsl_stepper_validation.txt' using (bin($3-$14,binwidth)):(1.0) smooth freq with boxes 

set table 'tmp_px'
plot 'gsl_stepper_validation.txt' using (bin($4-$15,binwidth)):(1.0) smooth freq with boxes 

set table 'tmp_py'
plot 'gsl_stepper_validation.txt' using (bin($5-$16,binwidth)):(1.0) smooth freq with boxes 

set table 'tmp_pz'
plot 'gsl_stepper_validation.txt' using (bin($6-$17,binwidth)):(1.0) smooth freq with boxes 

set table 'tmp_phi'
binwidth=0.0001
plot 'gsl_stepper_validation.txt' using (bin($9-$20,binwidth)):(1.0) smooth freq with boxes 

set table 'tmp_theta'
plot 'gsl_stepper_validation.txt' using (bin($10-$21,binwidth)):(1.0) smooth freq with boxes 

set table 'tmp_qop'
plot 'gsl_stepper_validation.txt' using (bin(($11-$22)/1000,binwidth)):(1.0) smooth freq with boxes 

unset table 
set logscale y
set title "Position resiudals"
set xlabel "Delta [mm]"
set output "position.png"
plot 'tmp_x' with boxes title "x", \
     'tmp_y' with boxes title "y", \
     'tmp_z' with boxes title "z", \

set title "Momentum resiudals"
set xlabel "Delta [MeV]"
set output "momentum.png"
plot 'tmp_px' with boxes title "px", \
     'tmp_py' with boxes title "py", \
     'tmp_pz' with boxes title "pz", \

set title "Direction resiudals"
set xlabel "Delta [rad]"
set xrange [-1e-3:1e-3]
set output "direction.png"
plot 'tmp_phi' with boxes title "phi", \
     'tmp_theta' with boxes title "theta", \

set title "q/p resiudals"
set xlabel "Delta [1/MeV]"
set output "qop.png"
plot 'tmp_qop' with boxes title "q/p", \
