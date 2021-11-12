function y= correlate1(wf1,wf2);                      % C. Bergevin 11.24.14
% correlate two 1-D real-valued row vectors (should work similar to Matlab's xcorr.m) 
% +++
error(nargchk(2, 2, nargin)), error(nargoutchk(0, 1, nargout))
if ~isvector(wf1) || ~isvector(wf2)
    error('Parameters must be vectors.')
end
% ensure they are row vectors
if (~isrow(wf1)),   wf1= wf1';  end
if (~isrow(wf2)),   wf2= wf2';  end
m = length(wf1); n = length(wf2);    % extract relevant dimensions

% create new arrays as needed for operation
g= wf2; % no flip for wf2 (unlike the convolution)
f= [zeros(1,n) wf1 zeros(1,n)];
NN= m+n-1;
for k=1:NN
    % Note: It took me awhile to get this code right!
    y(k)= sum(conj(f).*[zeros(1,k) g zeros(1,m-k+n)]); % shifted wf2
end
return