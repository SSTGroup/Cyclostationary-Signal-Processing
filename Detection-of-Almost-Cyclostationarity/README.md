# Description

Files implementing the multiple-hypothesis-test-procedure for the detection of almost-cyclostationarity described in the paper "Detection of almost-cyclostationarity: An approach based on a multiple hypothesis test" by Stefanie Horstmann, David Ramirez and Peter J. Schreier, Proc. Asilomar Conference on Signals, Systems and Computers, Pacific Grove, CA, USA, November 2017. 

Furthermore, the technique is extended by a filterbank implementation of the required resampling stage and a cycle period estimator as proposed in the paper "Joint Detection of Almost-Cyclostationary Signals and Estimation of their Cycle Period" by S. Horstmann, D. Ramirez, and P.J. Schreier, submitted in IEEE Signal Processing Letters, May 2017.

# Abstract
We propose a technique that jointly detects the presence of almost-cyclostationary (ACS) signals in wide-sense stationary (WSS) noise and provides an estimate of their cycle period.
Since the cycle period of an ACS process is not an integer, the approach is based on a combination of a resampling stage and a multiple hypothesis test, which deal separately with the fractional part and the integer part of the cycle period. The approach requires resampling the signal at many different rates, which is computationally expensive. For this reason we propose a filter bank structure that allows us to efficiently resample a signal at many different rates by identifying common interpolation stages among the set of resampling rates. 
# Contact
In case of questions, suggestions, problems etc. please send an email.

Stefanie Horstmann: <stefanie.horstmann@sst.upb.de>
