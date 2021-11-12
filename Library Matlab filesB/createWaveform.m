function [Out] = createWaveform(In);   % 2017.09.11   C. Bergevin

% [functional, though room for improvement/expansion]
% Use: Functional form of EXspecREP3.m to be used to generate canonical 1-D
% waveforms for various purposes (e.g., generate acoustic stimuli, test
% signal processing routines)

% Output: two columns w/ In.Npoints rows, 1st column is time and 2nd column
% is the signal waveform

% ------
% Stimulus Type Legend (i.e., In.stimT)
% stimT= 0 - non-quantized sinusoid
% stimT= 1 - quantized sinusoid
% stimT= 2 - one quantized sinusoid, one un-quantized sinusoid
% stimT= 3 - two quantized sinusoids
% stimT= 4 - click (i.e., an impulse)
% stimT= 5 - noise (uniform distribution)
% stimT= 6 - chirp (flat mag.)
% stimT= 7 - noise (Gaussian distribution; flat spectrum, random phase)
% stimT= 8 - exponentially decaying sinusoid (i.e., HO impulse response)
% stimT= 9 - decaying exponential (low-pass filter)

% % % -----------------------------------------------
% ---
% Set default values (if unspecified)
if (~isfield(In,'stimT')), In.stimT= 1; end   % see legend above {1}
if (~isfield(In,'SR')), In.SR= 44100; end   % sample rate used to collect waveform [Hz] {44100}
if (~isfield(In,'Npoints')), In.Npoints= 8192; end   % # of points for the waveform {8192}
% --- Params for stimT=1-3,8
if (~isfield(In,'f')), In.f= 2000; end  % Frequency (for waveforms w/ tones) [Hz] {2000}
if (~isfield(In,'ratio')), In.ratio= 1.22; end   % specify f2/f2 ratio (for waveforms w/ two tones) {1.22}
% --- Params for stimT=4 (noise)
if (~isfield(In,'CLKon')), In.CLKon= 1000; end   % index at which click turns 'on' (starts at 1)
if (~isfield(In,'CLKoff')), In.CLKoff= 1001; end   % index at which click turns 'off'
% --- Params for stimT=6 (chirp)
if (~isfield(In,'f1S')), In.f1S= 2000.0; end   % if a chirp (In.stimT=2) starting freq. [Hz] [freq. swept linearly w/ time]
if (~isfield(In,'f1E')), In.f1E= 4000.0; end  % ending freq.  (energy usually extends twice this far out)
if (~isfield(In,'cType')), In.cType= 0; end  % boolean re sweept type: 0-linear, 1-log {0}
% --- Params for stimT= 8
if (~isfield(In,'alpha')), In.alpha= 50; end  % decay const. for sinusoid {50?}
% --- Params for stimT=
% ** Note ** Other stimulus parameters (for stimT>=4) can be changed below
% % % -----------------------------------------------

% ---
dt= 1/In.SR;  % spacing of time steps
freq= [0:In.Npoints/2];    % create a freq. array (for FFT bin labeling)
freq= In.SR*freq./In.Npoints;
% ---
% quantize the freq. (so to have an integral # of cycles in time window)
df = In.SR/In.Npoints;
fQ= ceil(In.f/df)*df;   % quantized natural freq.
t=[0:1/In.SR:(In.Npoints-1)/In.SR];  % create an array of time points, Npoints long
% ---
% compute stimulus
if In.stimT==0 % non-quantized sinusoid
    signal= cos(2*pi*In.f*t);
    disp(sprintf(' \n *Stimulus* - (non-quantized) sinusoid, f = %g Hz \n', In.f));
    disp(sprintf('specified freq. = %g Hz', In.f));
elseif In.stimT==1     % quantized sinusoid
    signal= cos(2*pi*fQ*t);
    disp(sprintf(' \n *Stimulus* - quantized sinusoid, f = %g Hz \n', fQ));
    disp(sprintf('specified freq. = %g Hz', In.f));
    disp(sprintf('quantized freq. = %g Hz', fQ));
elseif In.stimT==2     % one quantized sinusoid, one un-quantized sinusoid
    signal= cos(2*pi*fQ*t) + cos(2*pi*In.ratio*fQ*t);
    disp(sprintf(' \n *Stimulus* - two sinusoids (one quantized, one not) \n'));
elseif In.stimT==3     % two quantized sinusoids
    fQ2= ceil(In.ratio*In.f/df)*df;
    signal= cos(2*pi*fQ*t) + cos(2*pi*fQ2*t);
    disp(sprintf(' \n *Stimulus* -  two sinusoids (both quantized) \n'));
elseif In.stimT==4     % click (kludge way to do it!)
    clktemp1= zeros(1,In.Npoints);
    clktemp2= ones(1,In.CLKoff-In.CLKon);
    signal= [clktemp1(1:In.CLKon-1) clktemp2 clktemp1(In.CLKoff:end)];
    disp(sprintf(' \n *Stimulus* - Click \n'));
elseif In.stimT==5     % noise (flat)
    signal= rand(1,In.Npoints); disp(sprintf(' \n *Stimulus* - Noise1 \n'));
elseif In.stimT==6     % chirp (flat)
    f1SQ= ceil(In.f1S/df)*df;      %quantize the start/end freqs. (necessary?)
    f1EQ= ceil(In.f1E/df)*df;
    % ---
    if (In.cType== 0),  fSWP= f1SQ + (f1EQ-f1SQ)*(In.SR/In.Npoints)*t;  % LINEAR sweep rate
    else % LOG sweep rate
        tauS= -t(end)/(log(f1EQ/f1SQ));
        fSWP= f1SQ* exp(-t/tauS);
    end
    signal = sin(2*pi*fSWP.*t); disp(sprintf(' \n *Stimulus* - Chirp \n'));
elseif In.stimT==7     % noise (Gaussian)
    Asize= In.Npoints/2 +1;
    % create array of complex numbers w/ random phase and unit magnitude
    for (n=1:Asize),    theta= rand*2*pi; N2(n)= exp(i*theta); end
    N2=N2';
    tNoise=irfft(N2);  % now take the inverse FFT
    % --- scale it down so #s are between -1 and 1 (i.e. normalize)
    if (abs(min(tNoise)) > max(tNoise)),    tNoise= tNoise/abs(min(tNoise));
    else    tNoise= tNoise/max(tNoise); end
    signal= tNoise';disp(sprintf(' \n *Noise* - Gaussian, flat-spectrum \n'));
elseif In.stimT==8 % exponentially decaying cos
    signal= exp(-In.alpha*t).*sin(2*pi*fQ*t); disp(sprintf(' \n *Exponentially decaying (quantized) sinusoid*  \n'));
elseif In.stimT==9 % decaying exponential
    signal= exp(-t*In.alpha); disp(sprintf(' \n *Decaying exponential*  \n'));
end

Out= [t; signal]';
return



