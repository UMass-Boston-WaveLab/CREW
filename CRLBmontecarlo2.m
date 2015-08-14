function [CRLBvsM, CRLBvsN, CRLBvsd, CRLBvsq, CRLBvsSNR] = CRLBmontecarlo2(reps)
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
M=31; %in general, make M be odd so CRLB calculation is simplified
N=5;
d=0.5;
SNR=20; %in dB: 10*log10(avg(|h|^2/|n|^2)) = 10*log10(N/sigma^2)
q=M+2/d; %2 wavelengths beyond sample region is default
mindk=pi/(8*M*d);

testM = 11:4:200;
testN = 1:20;
testd = 0.1:0.05:0.5;
testq = (M+1):(M+10/d);
testSNR = 1:30; %in dB
testdk = pi/(32*M*d):pi/(32*M*d):pi/(2*M*d);

temp=0;
CRLBvsM = zeros(size(testM));
CRLBvsN = zeros(size(testN));
CRLBvsd = zeros(size(testd));
CRLBvsq = zeros(size(testq));
CRLBvsSNR = zeros(size(testSNR));
CRLBvsdk = zeros(size(testdk));

rvsM = zeros(size(testM));
rvsN = zeros(size(testN));
rvsd = zeros(size(testd));
rvsq = zeros(size(testq));
rvsSNR = zeros(size(testSNR));
rvsdk = zeros(size(testdk));

for ii=1:length(testM)
    for jj = 1:reps
        [r, crlb]=chanestCRLB(N, testM(ii), d, q, 10^(SNR/10), mindk);
        temp=temp+crlb/reps;
        rvsM(ii) = rvsM(ii)+r/reps;
    end
    CRLBvsM(ii)=temp;
    temp=0;
    
end

for ii=1:length(testN)
    for jj = 1:reps
        [r, crlb]=chanestCRLB(testN(ii), M, d, q, 10^(SNR/10), mindk); %averaging
        temp=temp+crlb/reps;
        rvsN(ii) =rvsN(ii)+r/reps;
    end
    CRLBvsN(ii)=temp;
    temp=0;
    
end

for ii=1:length(testd)
    for jj = 1:reps
        [r, crlb]=chanestCRLB(N, M, testd(ii), q, 10^(SNR/10), mindk);
        temp=temp+crlb/reps;
         rvsd(ii) = rvsd(ii)+r/reps;
    end
    CRLBvsd(ii)=temp;
    temp=0;
   
end

for ii=1:length(testq)
    for jj = 1:reps
        [r, crlb]=chanestCRLB(N, M, d, testq(ii), 10^(SNR/10), mindk);
        temp=temp+crlb/reps;
        rvsq(ii) = rvsq(ii)+r/reps;
    end
    CRLBvsq(ii)=temp;
    temp=0;
    
end

for ii=1:length(testSNR)
    for jj = 1:reps
        [r, crlb]=chanestCRLB(N, M, d, q, 10^(testSNR(ii)/10), mindk);
        temp=temp+crlb/reps; 
        rvsSNR(ii) = rvsSNR(ii)+r/reps;
    end
    CRLBvsSNR(ii)=temp;
    temp=0;
    
end

for ii=1:length(testdk)
    for jj = 1:reps
        [r, crlb]=chanestCRLB(N, M, d, q, 10^(SNR/10), testdk(ii));
        temp=temp+crlb/reps; 
        rvsdk(ii) = rvsdk(ii)+r/reps;
    end
    CRLBvsdk(ii)=temp;
    temp=0;
    
end

%these are for debugging the default parameter choices
% figure;
% plot(testM, rvsM)
% xlabel('Measurement Length (wavelengths)')
% ylabel('Info Matrix Avg Rank')
% 
% figure;
% plot(testd, rvsd)
% xlabel('Space Between Samples')
% ylabel('Info Matrix Avg Rank')
% 
% figure;
% plot(testq, rvsq)
% xlabel('Number of Samples Predicted Ahead')
% ylabel('Info Matrix Avg Rank')


figure; 
plot(testM, real(CRLBvsM))
xlabel('Measurement Length (wavelengths)')
ylabel('Estimator Minimum Variance')
% ylim([0,1])

figure; 
plot(testN, real(CRLBvsN))
xlabel('Number of Scatterers')
ylabel('Estimator Minimum Variance')
% ylim([0,1])

figure; 
plot(testd, real(CRLBvsd))
xlabel('Space Between Samples (Wavelengths)')
ylabel('Estimator Minimum Variance')
% ylim([0,1])

figure; 
plot(testq-M, real(CRLBvsq))
xlabel(sprintf('Number of Samples Predicted Ahead (sample spacing = %.2f wavelengths)', d))
ylabel('Estimator Minimum Variance')
% ylim([0,1])

figure; 
plot(testSNR, real(CRLBvsSNR))
xlabel('Signal to Noise Ratio (dB)')
ylabel('Estimator Minimum Variance')
% ylim([0,1])

figure; 
plot(testdk, real(CRLBvsdk))
xlabel('Minimum Wavenumber Spacing (radians)')
ylabel('Estimator Minimum Variance')
% ylim([0,1])

end