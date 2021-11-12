function Out=makeKernel(In);
% create a kernel to be used for 2-D cross-correlation CB
% +++
% should be odd
if (~isfield(In,'resX')), In.resX= 25; end   % horizontal length of image (# of pixels) [400]
if (~isfield(In,'resY')), In.resY= 25; end   % vertical length of image (# of pixels) [300]
if (~isfield(In,'edgeT')), In.edgeT= 1; end   % hard (0) or soft (1) edge? [1]
if (~isfield(In,'Ampl')), In.Ampl= 255; end   % max. amplitude in image for each object [255] {requires In.edgeT= 1}
if (~isfield(In,'rollF')), In.rollF= 5; end   % baseline 'fall off' factor [5]
if (~isfield(In,'logS')), In.logS= 1; end   % compress to log scale (0-no,1-yes) [1]
% +++
Out.center= [ceil(In.resX/2) ceil(In.resY/2)];  % determine center point
% build up object image pixel-by-pixel (brute force!)
temp=[];
for j=1:In.resX
    for k=1:In.resY
        % center point
        if (j==Out.center(1) && k==Out.center(2))
            temp(j,k)= In.Ampl;
        % all other points
        else
            % edge falls off gradual (In.edgeT==1) or hard (In.edgeT==1)
            if In.edgeT==0
                temp(j,k)= 0;
            else
                dist= sqrt((j-Out.center(1))^2 + (k-Out.center(1))^2);
                temp(j,k)= In.Ampl/(In.rollF*dist);
            end
        end
    end
end
Out.array= temp';
% compress to log scale
if (In.logS==1),    Out.array= 10*log10(Out.array);  end
% plot image?
if (1==0),  imagesc(Out.array); colormap(gray); colorbar; caxis([0 max(max((Out.array)))]); end
return
