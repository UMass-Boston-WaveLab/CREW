function [ r, CRLB ] = chanestCRLB( N, M, d, q, SNR,mindk)
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

alpha = (1/sqrt(2))*(randn(1,N)+1i*randn(1,N)); %randn has variance 1, so 
                        %E(alpha^2) is 1 if alpha is constructed this way

dk=0;
tries=0;
while any(abs(dk)<mindk) && tries<10000
    theta = 2*pi*rand(1, N);
    psi = 2*pi*rand(1,N);
    k = 2*pi*(1+(1e-8)*cos(theta)).*cos(psi);
    dk = diff(sort(k));  %reject k vectors where some k are too close together
    tries=tries+1;
end
if any(dk<mindk)
    error(sprintf('Found no good k vector after %.0f tries, Md = %.2f, N= %.0f\n', tries, M*d, N))
end

%since SNR = N/sigma^2, sigma^2 = N/SNR.  But we normalize h to make
%|h|^2=1, so...
sigma = sqrt(1/SNR); %we never use sigma not-squared, but be consistent w/notation

%%create Hprime, the derivatives of h wrt parameters
Hprime = [0 exp(1i*k*q*d) 1i*exp(1i*k*q*d) 1i*alpha.*d*q.*exp(1i*k*q*d)];

%%create B, the CRLB for the parameters

%Xaa(i,j)=dh'/d(alpha_i)*dh/d(alpha_j)
Xaa = zeros(N);
for ii = 1:N
    for jj=1:N
        Xaa(ii,jj)=sum((exp(-j*k(ii)*d*(-(M-1)/2:(M-1)/2))).*( exp(j*k(jj)*d*(-(M-1)/2:(M-1)/2))));
    end
end

Xab = j*Xaa;
Xba = -j*Xaa;
Xbb = Xaa;

Xak = zeros(N);
for ii = 1:N
    for jj = 1:N
        Xak(ii,jj) = sum(alpha(jj)*j*d*(-(M-1)/2:(M-1)/2).*exp(j*(k(jj)-k(ii))*d*(-(M-1)/2:(M-1)/2)));
        if ii==jj 
            Xak(ii,jj)= 0;
        end
    end
end

Xbk = -j*Xak;

Xka = zeros(N);
for ii = 1:N
    for jj = 1:N
        Xka(ii,jj) = sum(conj(alpha(ii))*(-j)*d*(-(M-1)/2:(M-1)/2).*exp(j*(k(jj)-k(ii))*d*(-(M-1)/2:(M-1)/2)));
        if ii==jj
            Xka(ii,jj)= 0;
        end
    end
end

Xkb = -j*Xka;

Xkk = zeros(N);
for ii=1:N
    for jj=1:N
        Xkk(ii,jj)=sum(d^2*conj(alpha(ii))*alpha(jj)*(-(M-1)/2:(M-1)/2).^2.*exp(j*d*(-(M-1)/2:(M-1)/2)*(k(jj)-k(ii))));
    end
end



Binv = [M/sigma^4 zeros(1, 3*N); %should this be M^2?
        zeros(N,1) 2/(sigma^2)*real(Xaa) 2/(sigma^2)*real(Xab) 2/(sigma^2)*real(Xak);
        zeros(N,1) 2/(sigma^2)*real(Xba) 2/(sigma^2)*real(Xbb) 2/(sigma^2)*real(Xbk);
        zeros(N,1) 2/(sigma^2)*real(Xka) 2/(sigma^2)*real(Xkb) 2/(sigma^2)*real(Xkk)];
if cond(Binv)>10^10
    sprintf('oh noes');
end
r = rank(Binv);

    
CRLB=Hprime*(Binv\Hprime')/N^2;  

end

