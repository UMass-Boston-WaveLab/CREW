function [  ] = minvartest( N, M, d )
%test of minimum variance spectral estimation method
theta = rand(N,1);
k=2*pi*cos(theta);
h = sum(exp(j*kron(k,d*((1:M) - (M-1)/2))), 1);

powspec=minvar(h, floor(length(h)/3));

figure;
plot(powspec)


end

