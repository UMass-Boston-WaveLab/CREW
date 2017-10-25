function [] = MIvsCorr(maxR, deltar, M, N, L, bins, phi_evol, alpha_evol)
%performs a Monte Carlo simulation to demonstrate that mutual information
%can drop off more slowly than the correlation function, depending on what
%the underlying environmental parameters' spatial dependence is.
%maxR is the maximum distance between channel observations (in order to
%make sure there are enough pairs with each distance, r goes from -maxR to
%+maxR

%N is the number of scattering paths in a channel
%M is the number of channel realizations in a computation
%L is the characteristic length scale of the alphas' and phis' covariance
%functions
%bins = number of bins in histograms for prob dist of real and imaginary parts of h
%varp = variance of angle of arrival of a scattering path

r = -maxR:deltar:maxR; %in wavelengths
nr = length(r);
h = zeros(M, nr);

%% define covariance functions for alphas and phis
%this describes how fast channel changes
deltar = r;
% vara = 1/N; %normalizes total power of h to approx. 1
% K_alpha = exp(-deltar.^2/(2*L^2)); %"squared exponential" cov function - normalization happens below
% K_phi = exp(-deltar.^2/(2*L^2));

alphas = zeros(M,N,nr);
phis = zeros(M,N, nr);
%% make channel realizations
for ii = 1:M
% mus = 2*pi*rand(1,N); %randomly choose the mean for each phi
    parfor jj = 1:N
%        %alphas are complex gaussian process with mean 0, variance vara and cov function K_alpha
%        alphas(ii,jj,:) = sqrt(2/sqrt(pi))*sqrt(2*maxR/(nr*L))*ifft(fft(sqrt(1/2)*sqrt(vara)*(randn(size(r))+j*randn(size(r)))).*fft(K_alpha));
        alphas(ii,jj,:) = fill_alpha(alphas(ii,jj,:), alpha_evol);
        
%         %phis are real gaussian process with mean mus, variance varp, and
%         %cov function K_phi
%         phis(ii,jj,:) = sqrt(2/sqrt(pi))*sqrt(2*maxR/(nr*L))*ifft(fft(sqrt(varp)*randn(size(r))).*fft(K_phi));+mus(jj);  
        phis(ii,jj,:) = fill_phis(phis(ii,jj,:), phi_evol);
    end
     %make mth chan realization - alpha, phi, and r are funcs of r
    h(ii,:) = (1/N)*sum(squeeze(alphas(ii,:,:)).*exp(j*2*pi*sin(squeeze(phis(ii,:,:))).*repmat(r, N, 1)),1);
    %factor of 1/N makes average power be 1
    %estimated correlation func is corr(h)
    Rest = corr(h);
end



% %% make histograms vs. r
% 
% p_R = zeros(nr, nr, bins, bins);
% p_I = zeros(nr, nr, bins, bins);
% for ind1 = 1:nr
%     %loop through all the distances between points
%     parfor ind2 = 1:nr
%         %create p(re(h(r1)), re(h(r2))) by making a histogram of all the
%         %realizations' real and imaginary part values
%         p_R(ind1, ind2,:,:) = hist3(real([h(:,ind1) h(:,ind2)]), [bins, bins])/M;
%         p_I(ind1, ind2,:,:) = hist3(imag([h(:,ind1) h(:,ind2)]), [bins, bins])/M;
%         p_IR(ind1, ind2,:,:) = hist3([real(h(:,ind1)) imag(h(:,ind2))], [bins, bins])/M;
%         
%         %correlation func is ensemble average of h(x1)'*h(x2)
%         R(ind1,ind2) = (1/M)*h(:, ind1)'*h(:,ind2);
%     end
% end
% %if the channel process is wide-sense stationary vs. r, p_R and p_I
% %should be roughly constant vs. ind1 or ind2 but may vary wrt ind2-ind1
% 
% %% use probability distribution estimates to estimate mutual information
% I_R = zeros(nr, nr);
% I_I = zeros(nr, nr);
% I_I = zeros(nr, nr);
% for ii = 1:nr
%     for jj = 1:nr
%         I_R(ii,jj) = mutualinformation(squeeze(p_R(ii,jj,:,:)));
%         I_I(ii,jj) = mutualinformation(squeeze(p_I(ii,jj,:,:)));
%         I_IR(ii,jj) = mutualinformation(squeeze(p_IR(ii,jj,:,:)));
%     end
% end
% 


%% make histograms vs. r separation
if floor(nr/2)<nr/2
    ul = floor(nr/2)+1;
else
    ul = floor(nr/2);
end
upperR = zeros(ul, bins, bins);
upperI = zeros(ul, bins, bins);
upperIR = zeros(ul, bins, bins);

lowerR = zeros(floor(nr/2), bins, bins);
lowerI = zeros(floor(nr/2), bins, bins);
lowerIR = zeros(floor(nr/2), bins, bins);
for ii = 0:(ul-1) 
    upperR(ii+1,:,:) = hist3([reshape(real(h(:, 1:ul)),[],1) reshape(real(h(:, ii+(1:ul))),[],1)], [bins, bins]);
    upperI(ii+1,:,:) = hist3([reshape(imag(h(:, 1:ul)),[],1) reshape(imag(h(:, ii+(1:ul))),[],1)], [bins, bins]);
    upperIR(ii+1,:,:) = hist3([reshape(imag(h(:, 1:ul)),[],1) reshape(real(h(:, ii+(1:ul))),[],1)], [bins, bins]);
end
for ii = 1:floor(nr/2)
    lowerR(ii,:,:) = hist3([reshape(real(h(:, 1:floor(nr/2))),[],1) reshape(real(h(:, ((ul+1):nr)-ii)),[],1)], [bins, bins]);
    lowerI(ii,:,:) = hist3([reshape(imag(h(:, 1:floor(nr/2))),[],1) reshape(imag(h(:, ((ul+1):nr)-ii)),[],1)], [bins, bins]);
    lowerIR(ii,:,:) = hist3([reshape(imag(h(:, 1:floor(nr/2))),[],1) reshape(real(h(:, ((ul+1):nr)-ii)),[],1)], [bins, bins]);
end

p_deltaR = cat(1, lowerR, upperR)/(M*nr);
p_deltaI = cat(1, lowerI, upperI)/(M*nr);
p_deltaIR = cat(1, lowerIR, upperIR)/(M*nr);

for ii = 1:nr
    I_deltaR(ii) = mutualinformation(squeeze(p_deltaR(ii,:,:)), bins, M*nr);
    I_deltaI(ii) = mutualinformation(squeeze(p_deltaI(ii,:,:)), bins, M*nr);
    I_deltaIR(ii) = mutualinformation(squeeze(p_deltaIR(ii,:,:)), bins, M*nr);
end

figure;
plot(r,I_deltaR)
ylabel(sprintf('I(h(%.2f),h(x)), bits',r(ceil(nr/2))))
xlabel('\Delta_x (in \lambda)')
title(['Mutual Information of real(h(x)) and real(h(x+\Delta_x) vs. \Delta_x ' sprintf('char. length L=%.2f', L) ' \lambda'])
ylim([0 4])


% figure;
% plot(r, abs(R(ceil(nr/2),:)))
% ylabel(sprintf('R_h(%.2f,x)',r(ceil(nr/2))))
% xlabel('x (in \lambda)')
% title([sprintf('Correlation function (ensemble average) R_h(x_1, x_2) vs. x, char. length L=%.2f', L) ' \lambda'])


%real and imag parts are actually independent of each other
% figure;
% plot(r,I_IR(ceil(nr/2),:))
% ylabel(sprintf('I(h(%.2f),h(x)), bits',r(ceil(nr/2))))
% xlabel('x (in \lambda)')
% title([sprintf('Mutual Information of real(h(x)) and imag(h(%.2f', r(ceil(nr/2))) ' \lambda)) vs. x ' sprintf('char. length L=%.2f', L) ' \lambda'])
% 

figure;
plot(r, squeeze(abs(alphas(1,1,:))))
xlabel('x (in \lambda)')
title('Example of scattering path amplitude vs. position')

figure;
plot(r, squeeze(angle(alphas(1,1,:))))
xlabel('x (in \lambda)')
title('Example of scattering path initial phase vs. position')


figure;
plot(r, squeeze(phis(1,1,:)))
xlabel('x (in \lambda)')
title('Example of scattering path direction of arrival vs. position')
end

function [I] = mutualinformation(p, bins, samps)
%p is a 2D matrix that's the joint pdf of two random variables
%For a straightforward explanation: http://ai.stanford.edu/~gal/Research/Redundancy-Reduction/Neuron_suppl/node2.html
summand = p.*log2(p./(repmat(sum(p,1), size(p,1), 1).*repmat(sum(p,2),1, size(p,2))));
summand(isnan(summand)) = 0; %apply convention that 0*log(0) = 0
I = sum(sum(summand)) - bins/(2*samps*log2(2));
end

function [alphas] = fill_alpha(alphas, rho)
alphas(1) = (randn+j*randn)/sqrt(2);
for ii = 2:length(alphas)
    alphas(ii) = alphas(ii-1)+rho*(randn+j*randn)/sqrt(2); %complex
end
end

function [phis] = fill_phis(phis, rho)
phis(1) = 2*pi*rand;
for ii = 2:length(phis)
    phis(ii) = wrapToPi(phis(ii-1)+rho*randn);
end
end