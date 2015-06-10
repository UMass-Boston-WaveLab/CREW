function [CRLBvsM, CRLBvsN, CRLBvsd, CRLBvsq, CRLBvsSNR] = CRLBmontecarlo(reps, q)
%CRLBMONTECARLO Uses the Monte Carlo method to estimate the Cramer-Rao 
%Lower Bound on the variance of an unbiased estimator of the channel 
%transfer function at one location based on M spatial samples at other 
%locations.
%   Channel transfer function is assumed to be sum-of-sinusoids form
%   Not worried about estimating how many scatterers there are - that is a
%   separate problem
%   d (separation between samples) is specified in wavelengths
%   CRLB depends on k, alpha, which are random variables.  We can assume a
%   distribution for them and then calculate the actual expected value (this is
%   probably more general/correct) or we can Monte Carlo it.
%   We're going to do it the Monte Carlo way, since that's what Svantesson
%   and Swindlehurst do in their MIMO CRLB paper.
%   In this program we sweep through many values for various estimation 
%   parameters.

%default values
M=100;
N=10;
d=0.1;
SNR=20; %in dB: 10*log10(avg(|h|^2/|n|^2)) = 10*log10(N/sigma^2)

testM = 10:10:300;
testN = 1:30;
testd = 0.05:0.05:0.5;
testq = floor(0.25/d):floor(0.25/d):floor(10/d);
testSNR=13:30; %in dB

temp=0;
CRLBvsM = zeros(size(testM));
CRLBvsN = zeros(size(testN));
CRLBvsd = zeros(size(testd));
CRLBvsq = zeros(size(testq));
CRLBvsSNR = zeros(size(testSNR));

for ii=1:length(testN)
    for jj = 1:reps
        temp=temp+chanestCRLB(testN(ii), M, d, q, 10^(SNR/10))/reps; %averaging
    end
    CRLBvsN(ii)=temp;
    temp=0;
end

for ii=1:length(testM)
    for jj = 1:reps
        temp=temp+chanestCRLB(N, testM(ii), d, q, 10^(SNR/10))/reps;
    end
    CRLBvsM(ii)=temp;
    temp=0;
end

for ii=1:length(testd)
    for jj = 1:reps
        temp=temp+chanestCRLB(N, M, testd(ii), q, 10^(SNR/10))/reps;
    end
    CRLBvsd(ii)=temp;
    temp=0;
end

for ii=1:length(testq)
    for jj = 1:reps
        temp=temp+chanestCRLB(N, M, d, testq(ii), 10^(SNR/10))/reps;
    end
    CRLBvsq(ii)=temp;
    temp=0;
end

for ii=1:length(testSNR)
    for jj = 1:reps
        temp=temp+chanestCRLB(N, M, d, q, 10^(testSNR(ii)/10))/reps; 
    end
    CRLBvsSNR(ii)=temp;
    temp=0;
end
figure; 
plot(testM, real(CRLBvsM))
xlabel('Number of Measurements')
ylabel('Estimator Minimum Variance')

figure; 
plot(testN, real(CRLBvsN))
xlabel('Number of Scatterers')
ylabel('Estimator Minimum Variance')

figure; 
plot(testd, real(CRLBvsd))
xlabel('Space Between Samples (Wavelengths)')
ylabel('Estimator Minimum Variance')

figure; 
plot(testq, real(CRLBvsq))
xlabel('Number of Samples Predicted Ahead (sample spacing = lambda/10)')
ylabel('Estimator Minimum Variance')

figure; 
plot(testSNR, real(CRLBvsSNR))
xlabel('Signal to Noise Ratio (dB)')
ylabel('Estimator Minimum Variance')

end