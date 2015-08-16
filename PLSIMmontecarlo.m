function [PLSIMvsS, PLSIMvsN, PLSIMvsd, PLSIMvsP, PLSIMvsSNR,PLSIMvsf_d,PLSIMvsf_c] = PLSIMmontecarlo(reps)

%PLSIMMONTECARLO. This is built off of the CRLBMONTECARLO built by Prof.
%Kirby-Patel. I have used a fairly good base version of the channel for
%estimation as a default value set. If there is a really long run time then
%I will split each variable into it's own program and make it more
%parallel. -Eric 

%default values
L=40;
S    = 7;                 % # of Scatterers
d    = 0.3;               % Spacing between eavesdropper samples in wavelengths
N    = floor(L/d);        % # of sensor array samples
q    = 200;               % Number of samples ahead we attempt to predict
P    = 33;                % Number of complex sinusoids that make up the wireless channel
f_d  = 11000;             % doppler frequency
f_c  = 2400000;           % carrier frequency
SNR = 16;                   % Signal to Noise Ratio in dB.
averaging=100;              %effective SNR will be dB(averaging)+SNR
thresh= 100;                  %threshold for signal vs. noise space in the MUSIC algorithm

%test ranges
testS = 1:20;
testN = 100:20:400;
testd = 0.1:0.05:0.5;
testP = 5:2:50;
testSNR = 10:5:50; %in dB
%testf_d = 5000:4000:200000;
%testf_c = 1200000:5000:3600000;
testthresh=1:10:500;
testaveraging=1:10:1000;

%PL_Sim = PL_Security_Sim_pmusic();
q_maxS = zeros(size(testS));
q_maxN = zeros(size(testN));
q_maxd = zeros(size(testd));
% q_maxq = zeros(size(testq));
q_maxP = zeros(size(testP));
q_maxSNR = zeros(size(testSNR));
% q_maxf_c = zeros(size(testf_c));
% q_maxf_d = zeros(size(testf_d));
q_maxthresh=zeros(size(testthresh));
q_maxaveraging=zeros(size(testaveraging));
for kk = 1:reps   
    
    for ii=1:length(testS)
  
        [H, H_hat] = PL_Security_Sim_pmusic(testS(ii), N, d, q, P, SNR, thresh, averaging, 0); %tests the scenario
        err = abs((H-H_hat)/sqrt(mean(abs(H).^2))); % finds error percentage
        temp=find(err>0.05, 1);
        if isempty(temp)
            temp=length(H);
        end
        q_maxS(ii) = q_maxS(ii)+temp/reps;
    end
    
    for ii=1:length(testN)
        [H, H_hat] = PL_Security_Sim_pmusic(S, testN(ii), d, q, P,SNR, thresh, averaging,0); %tests the scenario
        err = abs((H-H_hat)/sqrt(mean(abs(H).^2))); % finds error percentage
        temp=find(err>0.05, 1);
        if isempty(temp)
            temp=length(H);
        end
        q_maxN(ii) = q_maxN(ii) + temp/reps;
    end
    
    for ii=1:length(testd)
        [H, H_hat] = PL_Security_Sim_pmusic(S, N, testd(ii), q, P,SNR, thresh, averaging,0); %tests the scenario
        
        err = abs((H-H_hat)/sqrt(mean(abs(H).^2))); % finds error percentage
        temp=find(err>0.05, 1);
        if isempty(temp)
            temp=length(H);
        end
        q_maxd(ii) = q_maxd(ii) + temp/reps;
    end
    
%    for ii=1:length(testq)
%         [H, H_hat] = PL_Sim(S, N, d, testq(ii), P,f_d,f_c,SNR); %tests the scenario
%         err = abs((H-H_hat)/sqrt(mean(abs(H).^2))); % finds error percentage
%         q_maxq(ii) = find(err>0.05, 'first');
%     end
    
    for ii=1:length(testP)
        [H, H_hat] = PL_Security_Sim_pmusic(S, N, d, q, testP(ii),SNR, thresh, averaging,0); %tests the scenario
        err = abs((H-H_hat)/sqrt(mean(abs(H).^2))); % finds error percentage
        temp=find(err>0.05, 1);
        if isempty(temp)
            temp=length(H);
        end
        q_maxP(ii) = q_maxP(ii) + temp/reps;
    end
    
%     for ii=1:length(testf_d)
%         [H, H_hat] = PL_Sim(testS(ii), N, d, q, P,testf_d(ii),f_c,SNR); %tests the scenario
%         err = abs((H-H_hat)/sqrt(mean(abs(H).^2))); % finds error percentage
%         q_maxf_d(ii) = find(err>0.05, 'first');
%     end
%     
%     for ii=1:length(testf_c)
%         [H, H_hat] = PL_Sim(testS(ii), N, d, q, P,f_d,testf_c(ii),SNR); %tests the scenario
%         err = abs((H-H_hat)/sqrt(mean(abs(H).^2))); % finds error percentage
%         q_maxf_c(ii) = find(err>0.05, 'first');
%     end
    for ii=1:length(testSNR)
        [H, H_hat] = PL_Security_Sim_pmusic(S, N, d, q, P,testSNR(ii), thresh, averaging,0); %tests the scenario
        err = abs((H-H_hat)/sqrt(mean(abs(H).^2))); % finds error percentage
        temp=find(err>0.05, 1);
        if isempty(temp)
            temp=length(H);
        end
        q_maxSNR(ii) = q_maxSNR(ii) + temp/reps;
    end
    for ii=1:length(testthresh)
        [H, H_hat] = PL_Security_Sim_pmusic(S, N, d, q, P,SNR, testthresh(ii), averaging,0); %tests the scenario
        err = abs((H-H_hat)/sqrt(mean(abs(H).^2))); % finds error percentage
        temp=find(err>0.05, 1);
        if isempty(temp)
            temp=length(H);
        end
        q_maxthresh(ii) = q_maxthresh(ii) + temp/reps;
    end
    for ii=1:length(testaveraging)
        [H, H_hat] = PL_Security_Sim_pmusic(S, N, d, q, P,SNR, thresh, testaveraging(ii),0); %tests the scenario
        err = abs((H-H_hat)/sqrt(mean(abs(H).^2))); % finds error percentage
        temp=find(err>0.05, 1);
        if isempty(temp)
            temp=length(H);
        end
        q_maxaveraging(ii) = q_maxaveraging(ii) + temp/reps;
    end
end
figure; 
plot(testS, q_maxS-N)
xlabel('Number of Scatterers')
ylabel('Average Max Prediction Past Array') %changed 'prediction' to 'agreement' because it starts at sample 1, not at end of measurement array

figure; 
plot(testN, q_maxN-testN)
xlabel('Number of Measurements')
ylabel('Average Max Prediction Length Past Array') %note that q_maxN has to have testN subtracted from it so that positive values indicate actual prediction (not estimation)

figure; 
plot(testd, q_maxd-N)
xlabel('Space Between Samples (Wavelengths)')
ylabel('Average Max Prediction Past Array')

figure; 
plot(testP, q_maxP-N)
xlabel('Number of component Sinusoids in rootMUSIC Estimate')
ylabel('Average Max Prediction Past Array')

% figure; 
% plot(testf_d, real(PLSIMvsf_d))
% xlabel(sprintf('Number of Samples Predicted Ahead (sample spacing = %.2f wavelengths)', d))
% ylabel('Estimator Minimum Variance')
% 
% figure; 
% plot(testf_c, real(PLSIMvsf_c))
% xlabel(sprintf('Number of Samples Predicted Ahead (sample spacing = %.2f wavelengths)', d))
% ylabel('Estimator Minimum Variance')

figure; 
plot(testSNR, q_maxSNR-N)
xlabel('Signal to Noise Ratio (dB)')
ylabel('Average Max Prediction Past Array')

figure; 
plot(testthresh, q_maxthresh-N)
xlabel('MUSIC Noise Subspace Threshold')
ylabel('Average Max Prediction Past Array')

figure; 
plot(testaveraging, q_maxaveraging-N)
xlabel('Number of Snapshots Averaged')
ylabel('Average Max Prediction Past Array')
