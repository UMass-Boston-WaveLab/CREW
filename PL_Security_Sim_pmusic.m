function [H, H_hat] = PL_Security_Sim_pmusic(S, N, d, q, P, SNR, smoo, Reps)
%%UMass Boston Physical Layer Security Channel Model
%Authors: Eric Brown, Clara Gamboa, Dr. K.C. Kerby-Patel
%
%       (SECTION 1) Variables
%       (SECTION 2) Parameters
%       (SECTION 3) Generated Signal
%       (SECTION 4) Noise
%       (SECTION 4) Specific Noise
%       (SECTION 5) PMusic Spectrum Analysis
%       (SECTION 6) Plotting Results
%       (SECTION 7) 
%--------------------------------------------


%%(SECTION 1)
%This is where we define the variables
 %S    = 7;                % # of Scatterers
 %N    = 200;               % # of sensor array samples
 %d    = 0.1;              % Spacing between eavesdropper samples in wavelengths
 %q    = 200;               % Number of samples ahead we attempt to predict
 Lamb = 1;                   % Wavelength = 1 (distances are normalized to the wavelength)
 t    = N+q;                 % is the total number of readings
 %P    = 60;                % Number of complex sinusoids that make up the wireless channel
 f_d  = 11000;             % doppler frequency
 f_c  = 2400000;           % carrier frequency
 %SNR = 13;                   % Signal to Noise Ratio.
 %Here we define a velocity vector for A

 
%---------------------------------------------


%% (SECTION 2) Parameters

    %In this section we define our wave characteristics. Phase angles and
    %thetas. In general A is used to describe the source node, B the receiving
    %node, and E the sensor array.


    pAS = 2*pi*rand(S, 1);      % pAS creates a vertical array of values representing
                                % the phase angles between S and source node (A)
    tAS = 2*pi*rand(S, 1);      % creates thetas for pAS

    pSE = 2*pi*rand(S, 1);      % pSE creates a vertical array of values representing
                                % the phase angles between S and sensor array (E)
    tSE = 2*pi*rand(S, 1);      % creates thetas for pSE
%---------------------------------------------


%% (SECTION 3)
% Generated signal for testing, generate what sensors N through N+q should
% be seeing.
% d is specified in terms of lambda, therefore, lambda = 1
% This defines the point on the wave as determined by distance from center
% of E. NOTE: This is the phase term not including any doppler effect

%AP is the Array Phase 
%ASC is the Array of Scatterers Continued???
%H is the full channel as seen by sensor array
k = (2*pi)/Lamb;

AP = repmat(pSE, [1,t]) + repmat(pAS, [1,t]) + kron(k*d*((1:t)-(N+1)/2),cos(tSE))+kron(k*(f_d/f_c)*d*((1:t)-(N+1)/2),(cos(tAS).*cos(tSE)));
ASC = exp(1i*AP);          
H = sum(ASC, 1);  

%---------------------------------------------

%% (SECTION 4) Noise 
% We use the pre-built matlab Gaussian White Noise
%Function to add AWGN to a given signal
Hn = zeros(size(H));
for hh = 1:Reps
    gWN = addgwn(H,SNR);
    Hn = (gWN/Reps)+Hn; 
end

%Hn = smoothing(Hn, smoo);

%---------------------------------------------

%% (SECTION 5)
% Here we use the pmusic spectrum analysis method. 
% k in the frequency and POW is the amplitude

% Variables
    x = Hn(1:N);        %sensor readings from our listening array
                        %p is the order of the linear preaditions(FIF filter)
                       %that predicts value of x
try
    [f,POW] = rootmusic(x,P);
catch
    f = zeros(P, 1);
    fprintf('root music has failed, P = %.0f S = %.0f  N = %.0f d= %.2f SNR = %.0f \n ', P, S, N, d, SNR ) 
end
% We have to make our sampling array into a matrix.(i.e 100 samples .1
% wavelength apart would result in a single row .1 - 10 intervals of .1, or
% 100 1 row columns
Samp = (1:t);

x_mat = repmat(Samp, size(f));

f_mat = repmat(f, size(Samp));

% Estimating complex amplitudes from the frequency
if any((abs(f)>0))
    z = exp(1i*f.');
    A = [];
    for ii = 1:N
        A = [A; z.^ii];
    end
    a = inv(A'*A)*A'*(Hn(1:N).');   
else    
    a = zeros(size(f));
end 
a_mat = repmat(a,size(Samp));
H_hat = sum(a_mat.*exp(1i*f_mat.*x_mat), 1);

%% (Section 6)
%Plotting the channel estimate vs the actual channel.
figure;
subplot(2,1,1)
plot(1:t-smoo,abs(H(1:t-smoo)),1:t-smoo,abs(H_hat(1:t-smoo)),1:t-smoo,abs(Hn(1:t-smoo)),'--'), grid
title 'Original Signal vs. rootMUSIC Estimate'
ylabel 'Amplitude (Linear)'
legend('Original signal','rootMUSIC Estimate', 'Signal with Noise')
subplot(2,1,2)
plot(1:t-smoo, 100*abs(H(1:t-smoo)-H_hat(1:t-smoo))/mean(abs(H(1:t-smoo)).^2)), grid
ylabel('Percent Error')
xlabel 'Sensors 1 through N+q' 
end