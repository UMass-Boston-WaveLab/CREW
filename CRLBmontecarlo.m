function [CRLBvsL, CRLBvsN, CRLBvsd, CRLBvsq, CRLBvsSNR] = CRLBmontecarlo(reps)
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
L=10;
N=10;
d=0.1;
SNR=13; %in dB: 10*log10(avg(|h|^2/|n|^2)) = 10*log10(N/sigma^2) unless h is normalized to 1 of course
q=floor(L/d)+2/d; %2 wavelengths beyond sample region is default

testL = 2:20;
testN = 1:30;
testd = 0.05:0.05:0.5;
testq = floor(L/d)+(floor(1/d):floor(10/d));
testSNR=7:30; %in dB

temp=0;
CRLBvsL = zeros(size(testL));
CRLBvsN = zeros(size(testN));
CRLBvsd = zeros(size(testd));
CRLBvsq = zeros(size(testq));
CRLBvsSNR = zeros(size(testSNR));

cnumvsL = zeros(size(testL));
cnumvsN = zeros(size(testN));
cnumvsd = zeros(size(testd));
cnumvsq = zeros(size(testq));
cnumvsSNR = zeros(size(testSNR));

for ii=1:length(testL)
    for jj = 1:reps
        [cnum, crlb]=chanestCRLB2(N, testL(ii), d, q, 10^(SNR/10));
        temp=temp+crlb/reps;
    end
    CRLBvsL(ii)=temp;
    temp=0;
    cnumvsL(ii) = max(cnum, cnumvsL(ii));
end

for ii=1:length(testN)
    for jj = 1:reps
        temp=temp+chanestCRLB2(testN(ii), L, d, q, 10^(SNR/10))/reps; %averaging
    end
    CRLBvsN(ii)=temp;
    temp=0;
    cnumvsN(ii) = max(cnum, cnumvsN(ii));
end

for ii=1:length(testd)
    for jj = 1:reps
        [cnum, crlb]=chanestCRLB2(N, L, testd(ii), q, 10^(SNR/10));
        temp=temp+crlb/reps;
    end
    CRLBvsd(ii)=temp;
    temp=0;
    cnumvsd(ii) = max(cnum, cnumvsd(ii));
end

for ii=1:length(testq)
    for jj = 1:reps
        [cnum, crlb]=chanestCRLB2(N, L, d, testq(ii), 10^(SNR/10));
        temp=temp+crlb/reps;
    end
    CRLBvsq(ii)=temp;
    temp=0;
    cnumvsq(ii) = max(cnum, cnumvsq(ii));
end

for ii=1:length(testSNR)
    for jj = 1:reps
        [cnum, crlb]=chanestCRLB2(N, L, d, q, 10^(testSNR(ii)/10));
        temp=temp+crlb/reps; 
        
    end
    CRLBvsSNR(ii)=temp;
    temp=0;
    cnumvsSNR(ii) = max(cnum, cnumvsSNR(ii));
end

figure;
plot(testL, 10*log10(cnumvsL))
xlabel('Measurement Length (wavelengths)')
ylabel('Info Matrix Max Condition Num [dB]')

figure;
plot(testd, 10*log10(cnumvsd))
xlabel('Space Between Samples')
ylabel('Info Matrix Max Condition Num [dB]')

figure;
plot(testq, 10*log10(cnumvsq))
xlabel('Number of Samples Predicted Ahead')
ylabel('Info Matrix Max Condition Num [dB]')


figure; 
plot(testL, real(CRLBvsL))
xlabel('Measurement Length (wavelengths)')
ylabel('Estimator Minimum Variance')
ylim([0,1])

figure; 
plot(testN, real(CRLBvsN))
xlabel('Number of Scatterers')
ylabel('Estimator Minimum Variance')
ylim([0,1])

figure; 
plot(testd, real(CRLBvsd))
xlabel('Space Between Samples (Wavelengths)')
ylabel('Estimator Minimum Variance')
ylim([0,1])

figure; 
plot(testq-floor(L/d), real(CRLBvsq))
xlabel(sprintf('Number of Samples Predicted Ahead (sample spacing = %.2f wavelengths)', d))
ylabel('Estimator Minimum Variance')
ylim([0,1])

figure; 
plot(testSNR, real(CRLBvsSNR))
xlabel('Signal to Noise Ratio (dB)')
ylabel('Estimator Minimum Variance')
ylim([0,1])

end