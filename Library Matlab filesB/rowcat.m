function B = rowcat (A, n)
% row-wise concatenates A with itself n times
% Matlab is such a crock o'shit that it doesn't
% have operators to do this.

% There must be a better way....
% B = A;
% for i = 1:(n-1)
%    B = [B A];
% end

% Now there is...
B = repmat (A,[1,n]);

