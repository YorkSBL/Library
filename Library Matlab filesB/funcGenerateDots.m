function Out=funcGenerateDots(In);  % 2014.05.14 CB

% [UPDATES IN PROGRESS/TESTING; still functional though]
% Purpose: Generates an array to be used as an image where there is a
% user-specified random object field (e.g., like a star field); 'intensity' of objects falls off like 1/r

% ---
% Input: structure 
% ---
% Output: structure w/ several fields
% ---
% Usage (e.g.,): >> x= funcGenerateDots; imagesc(x.array);
% ---
% Notes:
% o if In.AmplS= 0, all (max) intensities (need to test/verify)
% o if In.rollS= 0, all intensity falloff is identical (need to test/verify)
% if In.logS= 1, much "slower" fall-off (i.e., more diffuse)

% % ===============================================================================
% ---
%if ~(exist('In')),  clear; In=[];   end         % uncomment to run as a script
if (nargin < 1 | isempty (In)), In = [];           end;  % uncomment to run as a function
% +++
%if (~isstruct(In)), In= []; end   %
if (~isfield(In,'resX')), In.resX= 450; end   % horizontal length of image (# of pixels) [400]
if (~isfield(In,'resY')), In.resY= 300; end   % vertical length of image (# of pixels) [300]
if (~isfield(In,'Objects')), In.Objects= 12; end   % # of 'objects' to includes [12]
if (~isfield(In,'Ampl')), In.Ampl= 255; end   % max. amplitude in image for each object [255]
if (~isfield(In,'AmplS')), In.AmplS= 0.5; end   % scaling factor to turn of variable max. intensity [0.5?]
if (~isfield(In,'rollF')), In.rollF= 10; end   % baseline 'fall off' factor [5]
if (~isfield(In,'rollS')), In.rollS= 0.5; end   % scaling factor to turn of variable fall-off [0.5?]
if (~isfield(In,'logS')), In.logS= 1; end   % compress to log scale (0-no,1-yes) [1]
if (~isfield(In,'Noise')), In.Noise= 1; end   % add noise?
if (~isfield(In,'NoiseA')), In.NoiseA= 0.01; end   % noise amplitude scaling [0.01?]

% +++
Out.info= In;   % pass on parameters to output
Out.array= zeros(In.resY,In.resX);  % create base array

% loop to generate each object and (successively) add to base
for mm=1:In.Objects
    Out.coords(mm,:)= [ceil(rand(1,1)*In.resX) ceil(rand(1,1)*In.resY)]; % store coords.
    %Out.Ampl(mm,1)= rand* In.Ampl;  % randomized amplitude for given object
    %Out.roll(mm,1)= rand* In.rollF; % determine randomized falloff factor
    Out.Ampl(mm,1)= In.Ampl*(1+In.AmplS*rand);  % randomized amplitude for given object
    Out.roll(mm,1)= In.rollF*(1+In.rollS*rand); % determine randomized falloff factor
    
    % build up object image pixel-by-pixel (brute force!)
    temp=[];
    for j=1:In.resX
        for k=1:In.resY
            if (j==Out.coords(mm,1) && k==Out.coords(mm,2))
                temp(j,k)= Out.Ampl(mm,1);
            else
            dist= sqrt((j-Out.coords(mm,1))^2 + (k-Out.coords(mm,2))^2);
            temp(j,k)= Out.Ampl(mm,1)/(Out.roll(mm,1)*dist);
            end
        end
    end
    Out.array= Out.array+ temp';
end

% compress to log scale
if (In.logS==1),    Out.array= 10*log10(Out.array);  end

% add noise?
if (In.Noise==1), Out.array= Out.array+ In.NoiseA* In.Ampl* rand(In.resY,In.resX);  end

% plot image?
if (1==0),  imagesc(Out.array); colormap(gray); colorbar; caxis([0 max(max((Out.array)))]); end

return