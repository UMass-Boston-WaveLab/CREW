function [ CRLB ] = chanestCRLB( N, M, d, q)
%CHANESTCRLB Computes the Cramer-Rao Lower Bound on the variance of an
%unbiased estimator of the channel transfer function at one location based on
%M spatial samples.
%   Channel transfer function is assumed to be sum-of-sinusoids form
%   Not worried about estimating how many scatterers there are - that is a
%   separate problem
%   d (separation between samples) is specified in wavelengths
%   CRLB depends on k, alpha, which are random variables.  We can assume a
%   distribution for them and then calculate the actual expected value (this is
%   probably more general/correct) or we can Monte Carlo it.
%   We're going to do it the Monte Carlo way, since that's what Svantesson
%   and Swindlehurst do in their MIMO CRLB paper.




end

