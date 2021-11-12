function Out=findPeak2D(InB)    % 2014 CB
% Does a simple peak picking based upon scanning a WxW grid and checking to
% see if the center point is a max. (if so, store the coords.)

W= InB.W;  
%W=3;
Im= InB.im;
if (~isodd(W)), error('Specify odd #'); return; end   % req. odd # (so 'peak' is single point at center)
[M,N]= size(Im);    % determine array sizes

% loop through all possible overlap combinations for WxW search
indx=1;     % counter to keep track of getting a 'hit'
for ii=1:M-W-1
    for jj=1:N-W-1
        tp= Im(ii:ii+W-1,jj:jj+W-1);
        tpR= tp(:); % recast as a 1D vector
        % check if center point is a max
        flag= min(tp(ceil(W/2),ceil(W/2))>[tpR(1:floor(W^2/2));tpR(ceil(W^2/2)+1:W^2)]);
        if flag==1
            coords(indx,:)= [jj ii];
            indx= indx+ 1;
        end
    end
end

Out= coords;