function y= funcMakeCircle(P);    % 2017.01.22 CB (updated 2017.02.13)

% Function to make a two-tone "circle" (or annulus) "image", useful for
% things like Radon transforms and 2D FFT analysis
% ---
% Input params.
% P.N - # of "pixels" in the image (i.e., "image is NxN)
% P.r - radius of circle [pixels]
% P.fill - boolean to specify if the center should filled in (i.e., solid
% circle) or not (i.e., annulus)
% P.T - "thickness" of annulus [pixels]
% P.cent - center of circle [pixels]
% P.show - boolean to show resulting circle as an image
% P.invert - boolean to allow the "image" to be inverted
% P.noiseH - boolean to allow addition of Gaussian noise to "on" elements
% P.noiseB - boolean to allow addition of Gaussian noise to boundary locations
% P.offX - allows for controlled offset in x-position (re P.cent)
% P.offY - allows for controlled offset in y-position
% ---
% Notes
% o Code designed for pedagogical clarity, not computational
% efficiency/elegance. As such, there are two methods employed here: Method
% 1 is inefficient but spelled out, while Method 2 is more efficient but
% less intuitive (e.g., makes better use of Matlab's array nature). By
% default, Method 1 is commented out. Note also that as updates are made,
% they are only being made to Method 2 unless noted otherwise
% o Can turn into a script (rather than a function). Just swap the
% commmenting for the two conditionals below (and comment out the first
% line above if running as a script)
% o For Method 1, "noise" currently only added to intensity (not boundary)
% ---
%if ~(exist('P')),  clear; P=[];   end         % uncomment to run as a script
if (nargin < 1 | isempty (P)), P = [];           end;  % uncomment to run as a function
% ---
if (~isfield(P,'N')), P.N= 500; end  % # of pixels in (square) image length {500}
if (~isfield(P,'r')), P.r= P.N/4; end  % circle radius ["pixels"] {200}
if (~isfield(P,'fill')), P.fill= 0; end  % make a solid circle (i.e., not an annulus)? [boolean] {0}
if (~isfield(P,'T')), P.T= P.N/15; end  % circle thickness ["pixels"] {P.N/10} (Note: Only used if P.fill=0)
if (~isfield(P,'cent')), P.cent= floor(P.N/2); end  % circle center (make image in center by default)
if (~isfield(P,'show')), P.show= 1; end  % show resulting circle? [boolean] {0}
if (~isfield(P,'invert')), P.invert= 0; end  % boolean re inversion? {0}
if (~isfield(P,'offX')), P.offX= 0; end  % offset in x for center {0}
if (~isfield(P,'offY')), P.offY= 0; end  % offset in y for center {0}
% Note: For noise aspects below, degree of noise can be hacked further below in-line
if (~isfield(P,'noiseH')), P.noiseH= 0; end  % boolean re noise to "on" elements? {0}
if (~isfield(P,'noiseHs')), P.noiseHs= 0.25; end  % scaling factor re "on" element noise? {0.25} (reqs. P.noiseH=1)
if (~isfield(P,'noiseB')), P.noiseB= 0; end  % boolean re noise to boundary locations? {0}
if (~isfield(P,'noiseBs')), P.noiseBs= 0.05*P.r; end  % scaling factor re boundary location noise? {0.05*P.r} (reqs. P.noiseB=1)
% ---
% simple bookeeping checks
if (P.r>=P.N/2),  error('Specify a smaller radius');  end
if (P.r-P.T/2<=0),  disp('Thickness larger than circle!');  end
% % ===============================================================================
% % Method 1 (inefficient, but straightformward)
% % ---
% % create pixel array coords. (these are used to build up Z (the array
% % containing the circle)
% xx= [1:P.N];    yy= [1:P.N];
% % ---
% % create circle
% for nn=1:P.N
%     for mm=1:P.N
%         xT= xx(nn)-P.cent; yT= yy(mm)-P.cent; % shift "coords" to "center"
%         % conditional loop to specify (in binary fashion) is a given point
%         % is (or is not) to be set "on" (i.e., part of the circle)
%         if P.fill==1
%             % fill in the center
%             if (sqrt(xT^2+yT^2) <= P.r),    Z(nn,mm)=1;
%             else    Z(nn,mm)=0; end
%         else
%             % don't fill in the center (i.e., annulus)
%             if (sqrt(xT^2+yT^2) >= P.r-P.T/2) && (sqrt(xT^2+yT^2) <= P.r+P.T/2)
%                 Z(nn,mm)=1;
%             else    Z(nn,mm)=0; end
%         end
%     end
% end
% Z2= Z;
% ===============================================================================
% Method 2 (more efficient, but less straightformward)
[X,Y]= meshgrid([1:P.N]);
%X= X-P.cent; Y= Y-P.cent;    % shift "coords" to "center"
X= X-P.cent+P.offX; Y= Y-P.cent+P.offY;    % shift "coords" to "center" and offset if needed
Z2= sqrt(X.^2+Y.^2);        % make array of "radius" values
if (P.fill==1), Z2= Z2<=P.r+P.noiseB*P.noiseBs*randn(size(Z2));    % fill in the center
% don't fill in the center (i.e., annulus)
else    Z2= (Z2 >= P.r-P.T/2++P.noiseB*P.noiseBs*randn(size(Z2))) ...
        & (Z2 <= P.r+P.T/2++P.noiseB*P.noiseBs*randn(size(Z2))); end  
% ===============================================================================
if (P.noiseH==1)    % Add in intensity noise 
    Z2= double(Z2);     % need to convert from logical to double
    temp= find(Z2==1); Z2(temp)= Z2(temp)+ P.noiseHs*randn(numel(temp),1);   end
% ===============================================================================
% ---
if (P.invert==1),   Z2= ~(Z2==1);     end
% ---
if (P.show==1), figure(77); clf; imagesc(Z2); h1= colorbar; colormap bone; caxis([0 1]);
    xlabel('x [pixel #]'); ylabel('y [pixel #]'); h1.Label.String= 'Intensity'; end
% ---
y= Z2;
return