# Track Fitting
The track fitting algorithms estimate the track 
parameters. It is part of the pattern recognition/track 
reconstruction/tracking. We can run the 
track fitting algorithms, 
after we allocated all hits to single tracks with the 
help of a track finding algorithm. It is not necessary, 
that all points of a track are present.

Currently we have implementations for three different 
fitters:
* Kalman Filter
* GSF
* Global Chi-Square Fitter (GX2F) [wip]
?? is GSF a least square fit?
Eventhough all of them are least-squares fits, the 
concepts are quite different. Therefore, we should not 
expect identical results from all of them.


## Kalman Filter (KF)
The Kalman Filter is an iterative fitter. It 
successively combines measurments to obtain an 
estimate of the track parameters. The KF needs 
an estimate as a starting point. The procedure 
alternates between two methods:
1. Extrapolate the current state to the next 
layer.
2. Update the extrapolation using the measurment 
of the new layer.[^billoir]
The meaning of "this layer" and "the next layer" 
changes with the context. There are three 
different interpretations for this. The KF 
can give us those three interpretations as sets 
of track parameters:
* predicted: Uses "older" data to make the 
prediction.
* filtered: Uses the "current" data. Simply the 
weighted mean.
* smoothed: Uses the "future" data to predict 
the current parameters.



## GSF
(Benjamin)


## Global Chi-Square Fitter (GX2F) [wip]









[^billoir]: 
https://twiki.cern.ch/twiki/pub/LHCb/ParametrizedKalman/paramKalmanV01.pdf

