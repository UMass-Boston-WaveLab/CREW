function [PLSIMvsS, PLSIMvsN, PLSIMvsd, PLSIMvsP, PLSIMvsSNR,PLSIMvsf_d,PLSIMvsf_c] = PLSIMmontecarlo(reps)

%PLSIMMONTECARLO. This is built off of the CRLBMONTECARLO built by Prof.
%Kirby-Patel. I have used a fairly good base version of the channel for
%estimation as a default value set. If there is a really long run time then
%I will split each variable into it's own program and make it more
%parallel. -Eric 

%default values
S    = 6;                % # of Scatterers
N    = 200;               % # of sensor array samples
d    = 0.1;              % Spacing between eavesdropper samples in wavelengths
q    = 100;               % Number of samples ahead we attempt to predict
P    = 30;                % Number of complex sinusoids that make up the wireless channel
f_d  = 11000;             % doppler frequency
f_c  = 2400000;           % carrier frequency
SNR = 13;                   % Signal to Noise Ratio in dB.

%test ranges
testS = 1:20;
testN = 20:400;
testd = 0.05:0.05:0.5;
testP = 5:30;
testSNR = 5:30; %in dB
testf_d = 5000:4000:200000;
testf_c = 1200000:5000:3600000;

temp=0;

PLSIMvsS = zeros(size(testS));
PLSIMvsN = zeros(size(testN));
PLSIMvsd = zeros(size(testd));
PLSIMvsP = zeros(size(testP));
PLSIMvsSNR = zeros(size(testSNR));
PLSIMvsf_d = zeros(size(testf_d));
PLSIMvsf_c = zeros(size(testf_c));

PL_Sim = PL_Security_Sim_pmusic();
E_testS = [];
pred_Len_testS = [];
E_testN = [];
temp_length = []

for ii=1:length(testS)   
    PL_Sim(testS(ii), N, d, q, P,f_d,f_c,SNR); %tests the scenario 
    int jj = N;
    while err < 0.05 
        jj += 1;
        err = abs((H(jj)-H_hat(jj))/H(jj)); % finds error percentage
        
        % adds value of error to a vector. The length of this vector is the amount of samples out it can estimate
        E_testS(end+1)=err;
    end
   temp_length = length(E_testS);
   pred_Len_testS(end+1)=temp_length; 
    
end

for ii=1:length(testN)   
    PL_Sim(S, testN(ii), d, q, P,f_d,f_c,SNR); %tests the scenario 
    int jj = N;
    while err < 0.05 
        jj += 1;
        err = abs((H(jj)-H_hat(jj))/H(jj)); % finds error percentage
        
        % adds value of error to a vector. The length of this vector is the amount of samples out it can estimate
        E_testN(end+1)=err;
    end
   temp_length = length(E_testS);
   pred_Len_testS(end+1)=temp_length; 
    
end

for ii=1:length(testS)   
    PL_Sim(testS(ii), N, d, q, P,f_d,f_c,SNR); %tests the scenario 
    int jj = N;
    while err < 0.05 
        jj += 1;
        err = abs((H(jj)-H_hat(jj))/H(jj)); % finds error percentage
        
        % adds value of error to a vector. The length of this vector is the amount of samples out it can estimate
        E_testS(end+1)=err;
    end
   temp_length = length(E_testS);
   pred_Len_testS(end+1)=temp_length; 
    
end

for ii=1:length(testS)   
    PL_Sim(testS(ii), N, d, q, P,f_d,f_c,SNR); %tests the scenario 
    int jj = N;
    while err < 0.05 
        jj += 1;
        err = abs((H(jj)-H_hat(jj))/H(jj)); % finds error percentage
        
        % adds value of error to a vector. The length of this vector is the amount of samples out it can estimate
        E_testS(end+1)=err;
    end
   temp_length = length(E_testS);
   pred_Len_testS(end+1)=temp_length; 
    
end

for ii=1:length(testS)   
    PL_Sim(testS(ii), N, d, q, P,f_d,f_c,SNR); %tests the scenario 
    int jj = N;
    while err < 0.05 
        jj += 1;
        err = abs((H(jj)-H_hat(jj))/H(jj)); % finds error percentage
        
        % adds value of error to a vector. The length of this vector is the amount of samples out it can estimate
        E_testS(end+1)=err;
    end
   temp_length = length(E_testS);
   pred_Len_testS(end+1)=temp_length; 
    
end

for ii=1:length(testS)   
    PL_Sim(testS(ii), N, d, q, P,f_d,f_c,SNR); %tests the scenario 
    int jj = N;
    while err < 0.05 
        jj += 1;
        err = abs((H(jj)-H_hat(jj))/H(jj)); % finds error percentage
        
        % adds value of error to a vector. The length of this vector is the amount of samples out it can estimate
        E_testS(end+1)=err;
    end
   temp_length = length(E_testS);
   pred_Len_testS(end+1)=temp_length; 
    
end

for ii=1:length(testS)   
    PL_Sim(testS(ii), N, d, q, P,f_d,f_c,SNR); %tests the scenario 
    int jj = N;
    while err < 0.05 
        jj += 1;
        err = abs((H(jj)-H_hat(jj))/H(jj)); % finds error percentage
        
        % adds value of error to a vector. The length of this vector is the amount of samples out it can estimate
        E_testS(end+1)=err;
    end
   temp_length = length(E_testS);
   pred_Len_testS(end+1)=temp_length; 
    
end
figure; 
plot(testS, real(PLSIMvsS))
xlabel('Number of Measurements')
ylabel('Estimator Minimum Variance')

figure; 
plot(testN, real(PLSIMvsN))
xlabel('Number of Scatterers')
ylabel('Estimator Minimum Variance')

figure; 
plot(testd, real(PLSIMvsd))
xlabel('Space Between Samples (Wavelengths)')
ylabel('Estimator Minimum Variance')

figure; 
plot(testP, real(PLSIMvsP))
xlabel(sprintf('Number of Samples Predicted Ahead (sample spacing = %.2f wavelengths)', d))
ylabel('Estimator Minimum Variance')

figure; 
plot(testf_d, real(PLSIMvsf_d))
xlabel(sprintf('Number of Samples Predicted Ahead (sample spacing = %.2f wavelengths)', d))
ylabel('Estimator Minimum Variance')

figure; 
plot(testf_c, real(PLSIMvsf_c))
xlabel(sprintf('Number of Samples Predicted Ahead (sample spacing = %.2f wavelengths)', d))
ylabel('Estimator Minimum Variance')

figure; 
plot(testSNR, real(PLSIMvsSNR))
xlabel('Signal to Noise Ratio (dB)')
ylabel('Estimator Minimum Variance')

end