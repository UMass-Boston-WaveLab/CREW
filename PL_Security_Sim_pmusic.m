
%%UMass Boston Physical Layer Security Channel Model
%Authors: Eric Brown, Clara Gamboa, Dr. K.C. Kerby-Patel
%
%       (SECTION 1) Variables
%       (SECTION 2) Parameters
%       (SECTION 3) Generated Signal
%       (SECTION 4) Noise
%       (SECTION 5) PMusic Spectrum Analysis
%       (SECTION 6) 
%       (SECTION 7) 
%--------------------------------------------

%%  (SECTION 1)
%This is where we define the variables
 S = 10;                % # of Scatterers
 N = 100;               % # of sensor array samples
 d = 0.10;              % Spacing between eavesdropper samples in wavelengths
 q = 200;               % Number of samples ahead we attempt to predict
 Lamb = 1;              % Wavelength = 1 (distances are normalized to the wavelength)
 t = N+q;               % is the total number of readings
 P = 10;                % Number of complex sinusoids that make up the wireless channel
%---------------------------------------------

%% (SECTION 2)
%In this section we define our wave characteristics. Phase angles and
%thetas. In general A is used to describe the source node, B the receiving
%node, and E the sensor array.

pAS = 2*pi*rand(S, 1);  % pAS creates a vertical array of values representing
                        % the phase angles between S and source node (A)
tAS = 2*pi*rand(S, 1);  % creates thetas for pAS

pSE = 2*pi*rand(S, 1);  % pSE creates a vertical array of values representing
                        % the phase angles between S and sensor array (E)
tSE = 2*pi*rand(S, 1);  % creates thetas for pSE

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

% knum us the wavenumber
knum = kron(((2*pi)/Lamb)*d*((1:t)-(N+1)/2),cos(tSE));

% kd will be the doppler wave number, needs a doppler shift still
kd = kron(((2*pi)/Lamb)*d*((1:t)-(N+1)/2),cos(tAS)); 

AP = repmat(pSE, [1,t]) + repmat(pAS, [1,t]) + knum;
ASC = exp(1i*AP);
H = sum(ASC, 1);   
%---------------------------------------------

%% (SECTION 4)
% Here is where we define the noise added to the signal. 
% We have chosen to use a preset Matlab function, to generate our 
% Gaussian white noies. 
% Where gn is the function which adds white noise to our generated 
% signal H, using the specified signal to noise ratio defined 
% by the variable snr. 

% if more specific noise needs to be added later the function has the 
% ability to have the power of H in decibles relative to a watt defined
% as well as the imp is the specifies load impedance in ohms (default is 1)

%Hn = gn[ H ,;
%---------------------------------------------


%% (SECTION 5)
% Here we use the pmusic spectrum analysis method. 
% k in the frequency and POW is the amplitude

% Variables
    x = H(1:N);        %sensor readings from our listening array
    p = 10;            %p is the order of the linear preaditions(FIF filter)
                       %that predicts value of x

[k,POW] = rootmusic(x,P);

% We have to make our sampling array into a matrix.(i.e 100 samples .1
% wavelength apart would result in a single row .1 - 10 intervals of .1, or
% 100 1 row columns
Samp = (1:t);

x_mat = repmat(Samp, size(k));
%a_mat = repmat(POW, size(x));
k_mat = repmat(k, size(Samp));
% Estimating complex amplitudes from the frequency
z = exp(1i*k.');

A = [];
for ii = 1:N
    A = [A; z.^ii];
end



%a = inv(A)*H(1:S).;
a = inv(A'*A)*A'*(H(1:N).');
a_mat = repmat(a,size(Samp));


H_hat = sum(a_mat.*exp(1i*k_mat.*x_mat), 1);

%
plot(1:t,abs(H(1:t)),1:t,abs(H_hat),'--'), grid
title 'Original Signal vs. rootMUSIC Estimate'
xlabel 'Sensors 1 through N+q', ylabel 'Readings'
legend('Original signal','rootMUSIC Estimate')
