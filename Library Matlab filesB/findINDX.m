function y = findINDX(A,B)      % C. Bergevin 2.9.12

% feed in two arrays (A is nx1 and B is mx1 where m>=n) and the output will
% be the indicies of B whose elements are equal to those of A

% change dimensions if needbe
if ~iscolumn(A)
    A= A';
end

if ~iscolumn(B)
    B= B';
end


% make sure there is only a single column
if size(A,2)>1 || size(B,2)>1 || size(B,1)<size(A,1)
    fprintf('\n ** Error in array size ** \n');
    return
end


cnt= 1;
for nn=1:size(A,1)
    temp= find(A(nn)==B);
    if ~isempty(temp)
        indx(cnt)= temp;
        cnt= cnt+1;
    end
end

y= indx';

return