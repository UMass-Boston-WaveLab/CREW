function Px = minvar(x,p)

%MINVAR	Spectrum estimation using the minimum variance method.

%----

%USAGE	Px = minvar(x,p) 

%

%	The spectrum of a process x is estimated using the minimum

%	variance method (sometimes called the maximum likelihood

%	method).

%

%	x  :  Input sequence

%	p  :  Order of the minimum variance estimate - for short

%		sequences, p is typically about length(x)/3

%

%	The spectrum estimate is returned in Px using a dB scale.

%

%---------------------------------------------------------------

% copyright 1996, by M.H. Hayes.  For use with the book 

% "Statistical Digital Signal Processing and Modeling"

% (John Wiley & Sons, 1996).

%---------------------------------------------------------------

%

   x = x(:);

   R = covar(x,p);

   [v,d]=eig(R);

   U = diag(inv(abs(d)+eps));

   V  = abs(fft(v,1024)).^2;

   Px = 10*log10(p)-10*log10(V*U);

end