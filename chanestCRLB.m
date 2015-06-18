function [ CRLB ] = chanestCRLB( N, M, d, q, SNR)
%CHANESTCRLB Computes the Cramer-Rao Lower Bound on the variance of an
%unbiased estimator of the channel transfer function at one location based on
%M spatial samples.
%   Channel transfer function is assumed to be sum-of-sinusoids form
%   Not worried about estimating how many scatterers there are - that is a
%   separate problem
%   d (separation between samples) is specified in wavelengths
%   CRLB depends on k, alpha, which are random variables.  This computes
%   the CRLB for one particular channel realization - it is called by a
%   wrapper program that runs a Monte Carlo simulation to get an average
%   value.

%%generate random values that we can use to construct the channel: alpha_n,
%theta_n, psi_n, k_n
theta = 2*pi*rand(1, N);
psi = 2*pi*rand(1,N);
alpha = (1/sqrt(2))*(randn(1,N)+1i*randn(1,N)); %randn has variance 1, so 
                        %E(alpha^2) is 1 if alpha is constructed this way

k = 2*pi*(1+(1e-8)*cos(theta)).*cos(psi);

%since SNR = N/sigma^2, sigma^2 = N/SNR.  But we normalize h to make
%|h|^2=1, so...
sigma = sqrt(1/SNR); %we never use sigma not-squared, but be consistent w/notation

%%create Hprime, the derivatives of h wrt parameters
Hprime = 1/(sqrt(N))*[0 exp(1i*k*q*d) 1i*exp(1i*k*q*d) 1i*alpha.*d*q.*exp(1i*k*q*d)];

%%create B, the CRLB for the parameters
%Binv is made of different combinations of DR, DI, and Dk except for
%Binv(1,1).
DR = zeros(M,N);
DRrow = exp(1i*k*d); %k is a row vector
for ii = 1:M
    DR(ii,:) = 1/(sqrt(N))*DRrow.^ii;
end

DI = 1i*DR;

Dk = zeros(M,N);
for ii=1:M
    Dk(ii,:) = (1/sqrt(N))*1i*ii*d*alpha.*DRrow.^ii;
end

Binv = [M/sigma^4 zeros(1, 3*N);
        zeros(N,1) 2/(sigma^2)*real(DR'*DR) 2/(sigma^2)*real(DR'*DI) 2/(sigma^2)*real(DR'*Dk);
        zeros(N,1) 2/(sigma^2)*real(DI'*DR) 2/(sigma^2)*real(DI'*DI) 2/(sigma^2)*real(DI'*Dk);
        zeros(N,1) 2/(sigma^2)*real(Dk'*DR) 2/(sigma^2)*real(Dk'*DI) 2/(sigma^2)*real(Dk'*Dk)];
if cond(Binv)>10^5
    sprintf('oh noes');
end
    
    
    
CRLB=Hprime*(Binv\Hprime');
    
end

