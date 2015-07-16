function [ CRLB ] = chanestCRLB( N, L, d, q, SNR)
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
M=floor(L/d);
DRR = zeros(N,N);
for ii = 1:N
    for jj=1:N
        if ii==jj
            DRR(ii,jj) = M;
        else
            DRR(ii,jj) = dot(exp(j*k(ii)*((1:M)-(M-1)/2)*d), exp(-j*k(jj)*((1:M)-(M-1)/2)*d));
        end
    end
end

DRI = 1i*DRR;
DIR = conj(DRI);
DRk=zeros(size(DRR));
DIk = DRk;
DkR = DRk;
DkI = DRk;

Dkk = zeros(N,N);
for ii=1:N
    for jj = 1:N
    Dkk(ii,jj) = dot(d*(((1:M)-(M-1)/2).^2)*conj(alpha(ii)).*exp(j*k(ii)*((1:M)-(M-1)/2)*d),...
        d*(((1:M)-(M-1)/2).^2)*alpha(jj).*exp(-j*k(jj)*((1:M)-(M-1)/2)*d));
    end
end

Binv = [M/sigma^4 zeros(1, 3*N);
        zeros(N,1) 2/(sigma^2)*real(DRR) 2/(sigma^2)*real(DRI) DRk;
        zeros(N,1) 2/(sigma^2)*real(DIR) 2/(sigma^2)*real(DRR) DIk;
        zeros(N,1) DkR DkI 2/(sigma^2)*real(Dkk)];
if cond(Binv)>10^10
    sprintf('oh noes');
end
    
    
    
CRLB=Hprime*(Binv\Hprime');  %having a lot of trouble with B being "nearly singular" here
    
end

