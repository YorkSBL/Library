function Out=crosscorrIm(Ke,Im);        % 2014 C. Bergevin
% Does a simple cross-correlation analysis for two 2-D arrays, the 'kernel'
% (Ke) and the 'image' (Im); the dimensions of Ke should be less than Im

% determine array sizes
[m,n]= size(Ke);
[M,N]= size(Im);

% pad Im with zeroes on edges (so to allow for all possible overlap combinations; 'borrowed' older padding routine from Image Processing toolbox)
ImP= padarray(Im,[m-1 n-1],0,'both');

% loop through all possible overlap combinations
for ii=1:M+(m-1)
    for jj=1:N+(n-1)
        % determine appropriate sub-array from Im 
        corrVal(ii,jj)= sum(sum(Ke.*ImP(ii:ii+m-1,jj:jj+n-1)));
    end
end

Out= corrVal;