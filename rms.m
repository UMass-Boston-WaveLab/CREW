function [ rms ] = rms( varargin )
%Computes the root mean square value of the first argument, along the
%dimension specified by the second argument.

if length(varargin)==1
    X=varargin{1};
    rms = sqrt(sum(X.^2)/length(X));
elseif length(varargin)==2
    X=varargin{1};
    dim=varargin{2};
    if size(X,dim)~=0
        rms=sqrt(sum(X.^2,dim)/size(X,dim));
    else
        error('X has 0 length along specified dimension');
    end
else
    error('Wrong number of input arguments to rms()');
end


end

