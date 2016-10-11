set ylabel "Tracks"
set grid
set terminal pngcairo size 500,500 enhanced font 'Verdana,10'

unset logscale y
binwidth=0.001
bin(x,width)=width*floor(x/width)

set table 'tmp_x'
plot 'gsl_stepper_validation.txt' using (bin($1-$17,binwidth)):(1.0) smooth freq with boxes 

set table 'tmp_y'
plot 'gsl_stepper_validation.txt' using (bin($2-$18,binwidth)):(1.0) smooth freq with boxes 

set table 'tmp_z'
plot 'gsl_stepper_validation.txt' using (bin($3-$19,binwidth)):(1.0) smooth freq with boxes 

set table 'tmp_px'
plot 'gsl_stepper_validation.txt' using (bin($4-$20,binwidth)):(1.0) smooth freq with boxes 

set table 'tmp_py'
plot 'gsl_stepper_validation.txt' using (bin($5-$21,binwidth)):(1.0) smooth freq with boxes 

set table 'tmp_pz'
plot 'gsl_stepper_validation.txt' using (bin($6-$22,binwidth)):(1.0) smooth freq with boxes 

set table 'tmp_phi'
binwidth=0.0001
plot 'gsl_stepper_validation.txt' using (bin($9-$25,binwidth)):(1.0) smooth freq with boxes 

set table 'tmp_theta'
plot 'gsl_stepper_validation.txt' using (bin($10-$26,binwidth)):(1.0) smooth freq with boxes 

set table 'tmp_qop'
plot 'gsl_stepper_validation.txt' using (bin(($11-$27)/1000,binwidth)):(1.0) smooth freq with boxes 

set table 'rel_x'
plot 'gsl_stepper_validation.txt' using (bin(($1-$17)/$17,binwidth)):(1.0) smooth freq with boxes 

set table 'rel_y'
plot 'gsl_stepper_validation.txt' using (bin(($2-$18)/$18,binwidth)):(1.0) smooth freq with boxes 

set table 'rel_z'
plot 'gsl_stepper_validation.txt' using (bin(($3-$19)/$19,binwidth)):(1.0) smooth freq with boxes 

set table 'rel_px'
plot 'gsl_stepper_validation.txt' using (bin(($4-$20)/$20,binwidth)):(1.0) smooth freq with boxes 

set table 'rel_py'
plot 'gsl_stepper_validation.txt' using (bin(($5-$21)/$21,binwidth)):(1.0) smooth freq with boxes 

set table 'rel_pz'
plot 'gsl_stepper_validation.txt' using (bin(($5-$22)/$22,binwidth)):(1.0) smooth freq with boxes 

set table 'rel_phi'
plot 'gsl_stepper_validation.txt' using (bin(($9-$25)/$25,binwidth)):(1.0) smooth freq with boxes 

set table 'rel_theta'
plot 'gsl_stepper_validation.txt' using (bin(($10-$26)/$26,binwidth)):(1.0) smooth freq with boxes 

set table 'rel_qop'
plot 'gsl_stepper_validation.txt' using (bin(($11-$27)/$27,binwidth)):(1.0) smooth freq with boxes 

set table 'cov_phi'
plot 'gsl_stepper_validation.txt' using (bin(($14-$30)/$30,binwidth)):(1.0) smooth freq with boxes 

set table 'cov_theta'
plot 'gsl_stepper_validation.txt' using (bin(($15-$31)/$30,binwidth)):(1.0) smooth freq with boxes 

unset table 
set logscale y
set title "Position resiudals"
set xlabel "Delta [mm]"
#set xrange [-10:10]
set output "position.png"
plot 'tmp_x' with boxes title "x", \
     'tmp_y' with boxes title "y", \
     'tmp_z' with boxes title "z", \

set title "Momentum resiudals"
set xlabel "Delta [MeV]"
#set xrange [-10:10]
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

set title "relative differences"
set xlabel "delta/reference"
set output "relative.png"
plot 'rel_x' with boxes title "x", \
	 'rel_y' with boxes title "y", \
	 'rel_z' with boxes title "z", \
	 'rel_px' with boxes title "px", \
	 'rel_py' with boxes title "py", \
	 'rel_pz' with boxes title "pz", \
	 'rel_phi' with boxes title "phi", \
	 'rel_theta' with boxes title "theta", \
	 'rel_qop' with boxes title "q/p", \

set title "relative differences"
set xlabel "delta/reference"
set output "covariance.png"
plot 'cov_phi' with boxes title "cov(phi)", \
	 'cov_theta' with boxes title "cov(theta)", \

	 