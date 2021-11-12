function b = randINDX(x)
% randomize input (Nx1) array
%ix = randperm(n);
%b = A(ix);
b= x(randperm(numel(x)));
return