function [ Hn ] = gn( H , snr )
%gn - Guassian White Noise function 
 % for our noise we are using the Matlab function y = awgn (x,snr) 
 % awgn(x,snr) adds white Gaussian noise to the vector signal x.
 % the scalar snr specifies the singal to noise ratio per 
 % sample, in dB. If x is complex, awgn adds complex noise.

 
                                %wgm function generates an m by n matrix of white Gausssan noise
%hAverage = mean (H.');
%yTemp = wgn(m,n,p);            % y = wgn(m,n,p) generates an m-by-n matrix of white Gaussian noise. 
%nAverage = mean (yTemp.');
hAverage = mean(H.')
L=100000; %Sample length for the random signal
mu=0;
sigma=2;
X=sigma*randn(L,1)+mu;
snr = 
Hn = awgn(x,snr);     
 
end

