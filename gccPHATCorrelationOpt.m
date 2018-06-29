function [c,lags] = gccPHATCorrelationOpt(a,b)
% This function calculates the unscaled cross-correlation of 2
% vectors of the same length. The output length(c) is
% length(a)+length(b)-1. It is a simplified function of xcorr
% function in matlabR12 using the definition:
% c(m) = E[a(n+m)*conj(b(n))] = E[a(n)*conj(b(n-m))]
%
a = a(:); % convert a to column vector
b = b(:); % convert b to column vector

if length(a) > length(b)
    b = [b;zeros(length(a) - length(b),1)];
else
    a = [a;zeros(length(b) - length(a),1)];
end

M = length(a); % same as length(b)
maxlag = M-1; % maximum value of lag
lags = (-maxlag:maxlag)'; % vector of lags
A = fft(a,2^nextpow2(2*M-1)); % fft of A
A = A./abs(A);
B = fft(b,2^nextpow2(2*M-1)); % fft of B
B = B./abs(B);
c = ifft(A.*conj(B)); % crosscorrelation
%
% Move negative lags before positive lags.
%
c = [c(end-maxlag+1:end,1);c(1:maxlag+1,1)];
%
% Return row vector if a, b are row vectors.
%
[nr nc]=size(a);
if(nr>nc)
c=c.';
lags=lags.';
end
% End of function file.