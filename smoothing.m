function [ xsmooth ] = smoothing( x, w )
%implements smoothing

for ii = 1:(length(x)-w)
    xsmooth(ii)=mean(x(ii:(ii+w)));
end


end

