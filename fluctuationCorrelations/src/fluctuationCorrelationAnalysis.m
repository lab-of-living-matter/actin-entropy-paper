function [nfc,dfc] = fluctuationCorrelationAnalysis(stackPath)
if exist('stackPath','var') ~= 1 || isempty(stackPath)
    stackPath = [];
end

[expt,alignspec,stack] = defineParamsFlucCorrManuscript(stackPath);
[FFTAlignmentData] = FFTAlignmentStack_Manuscript(stack,alignspec);
[nfc] = directorFluctuationCorrelations(expt,FFTAlignmentData);
[dfc] = densityFluctuationCorrelations(expt,alignspec,stack,FFTAlignmentData);
end

function [expt,alignspec,stack] = defineParamsFlucCorrManuscript(stackPath)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script defines all of the relevent experimental and analysis
% parameters.
%                      
% Parameters specified:
%   
%       stack       = read tiff stack
%
%   Experimental - Tiff stack information & metadata         
%       expt          = structure containing all of the experimental
%                      parameters
%       expt.path     = tiff stack local path
%       expt.filename = tiff stack file name
%       expt.dscl     = image distance scale [=] µm/pixel
%       expt.tscl     = image time scale [=] seconds/frame
%       expt.timeind  = tiff stack frame index
%       expt.time     = tiff stack time vector (time index * time scale)
%       expt.lastim   = last image index to be analyzed (0=all images)
%       expt.stackdim = tiff stack dimensions [N x M x P]
%
%       expt.cp       = user defined center point for analyzing radially
%                      symmetric phenomena (set to [nan,nan] otherwise)                 
%       expt.maxr     = maximum radius (µm) for analyzing radially
%                      symmetric phenomena (set to [nan] otherwise) 
%
%   FFT Alignment
%       alignspec             = structure containing all of the FFT
%                               Alignment analysis parameters   
%       alignspec.winsize     = (pivspec.winsize)+1; %default
%       alignspec.overlap     = ((pivspec.winsize)*pivspec.overlap)/alignspec.winsize;
%       alignspec.st          = 2; % default
%       alignspec.checkpoint  = 0; %default
%       alignspec.mask_method = 1; %Global (val = 1) uses a circle for every window. Local (val = 2) uses a local threshold for each subwindow. (Default = 1)
%       alignspec.figures     = 1; 
%
%       winsize     = Size of the sub windows. Must be odd. (Default is 33)
%       overlap     = How much the windows overlap. Value should be between 0
%                     (complete overlap) and 1 (no overlap). (Default is
%                     0.5)
%       st          = order parameter spacing. Considers a window of size
%                     2*st+1 to compare vectors to. Depending on the
%                     overlap size, you'll want to change this. (Default is
%                     2)
%       checkpoint  = threshold value above which to do the caculation. The
%                     sum of the image intensity in the window needs to be
%                     above this value for the routine to calculate a vector
%                     for the given window (Default is 0 - i.e. every
%                     window)
%       mask_method = Determines the method to mask the FFT when
%                     determining the moment calculation. Global (val = 1) uses
%                     a circle for every window. Local (val = 2) uses a
%                     local threshold for each subwindow. (Default = 1)
%       figures     = Plot figures. 0 = No; 1 = Yes; (Default = 1)
%
%
%
% All rights and permissions belong to:
% Ian Linsmeier
% ian.linsmeier@yale.edu
% Lab of Living Matter (Murrell Lab)
% 4/15/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Experiment Parameters
if exist('stackPath')~=1 || isempty(stackPath)
    [FileName,PathName] = uigetfile('*.tif');
    stack = stackread([PathName,filesep,FileName]);
else
    [PathName,FileName] = fileparts(fullfile(stackPath));
    stack = stackread([PathName,filesep,FileName]);
end
stackdim = size(stack);
if length(stackdim) == 2, stackdim(3) = 1; end

expt.path = PathName;
expt.filename = FileName;
expt.dscl = 0.102;
expt.tscl = 1;
expt.timeind = linspace(0,stackdim(3)-1,stackdim(3));
expt.time = expt.timeind.*expt.tscl;
expt.lastim = 0;
expt.stackdim = stackdim;
expt.cp = [nan,nan];
expt.maxr = (min([stackdim(1),stackdim(2)])-200).*expt.dscl; %um


%% FFT Alignment Parameters
winsize = 33; % pixels
winspc = (winsize-1)/2; % ~50% overlap
alignspec.winsize = winsize;
alignspec.overlap = 1-(winspc/winsize);
alignspec.st = 2; % default
alignspec.checkpoint = 0; %default
alignspec.mask_method = 1; %Global (val = 1) uses a circle for every window. Local (val = 2) uses a local threshold for each subwindow. (Default = 1) 
alignspec.figures = 0;

end

function [FFTAlignmentData] = FFTAlignmentStack_Manuscript(stack,alignspec)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This routine produces a vectorfield of alignment directions using small
% subwindows in the real space image.  Images should be grayscale images.
% Angles determined are oriented as follows:
%
%                            ^  180?
%                            |
%                            |
%                            ------>  90?
%                            |
%                            |
%                            v  0?
%
% Inputs:
%       im          = the image to analyze. Can be any bit-depth.
%       winsize     = Size of the sub windows. Must be odd. (Default is 33)
%       overlap     = How much the windows overlap. Value should be between 0
%                     (complete overlap) and 1 (no overlap). (Default is
%                     0.5)
%       st          = order parameter spacing. Considers a window of size
%                     2*st+1 to compare vectors to. Depending on the
%                     overlap size, you'll want to change this. (Default is
%                     2)
%       checkpoint  = threshold value above which to do the caculation. The
%                     sum of the image intensity in the window needs to be
%                     above this value for the routine to calculate a vector
%                     for the given window (Default is 0 - i.e. every
%                     window)
%       mask_method = Determines the method to mask the FFT when
%                     determining the moment calculation. Global (val = 1) uses
%                     a circle for every window. Local (val = 2) uses a
%                     local threshold for each subwindow. (Default = 1)
%       figures     = Plot figures. 0 = No; 1 = Yes; (Default = 1)
%
%
% Outputs:
%       Everything is stored in a structure file FFTAlignmentData
%       anglemat     = A matrix of the calculated alignment direction in degrees. 
%       ordermat     = A matrix of the calculated order parameter;
%       pos          = row and column positions for each vector
%       vec          = vector direction components [coldisp rowdisp]
%       hist_bins    = bin values for the histogram of angles
%       hist_val     = # of angles in each bin
%       his_val_norm = normalized values of the histogram
%       parameters   = parameters used in the routine
%
%
% All rights and permissions belong to:
% Patrick Oakes
% poakes@uchicago.edu
% Gardel Laboratory
% 01/10/2013
%
% Updated by: Ian Linsmeier
% Updated last: 4/20/17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Select image stack to be analyzed
if exist('stack')~=1 || isempty(stack)
    [FileName,PathName] = uigetfile('*.tif');
    stack = stackread([PathName,'\',FileName]);
    stackdim = size(stack);
    if length(stackdim) == 2, stackdim(3) = 1; end
else
    stackdim = size(stack);
    if length(stackdim) == 2, stackdim(3) = 1; end
end

%% Run FFTAlignment on each image frame
for i=1:stackdim(3)
    FFTAlignmentData(i) = FFTAlignment(stack(:,:,i),alignspec.winsize,alignspec.overlap,alignspec.st,alignspec.checkpoint,alignspec.mask_method,alignspec.figures);
    disp(['Finished Frame ',num2str(i)]);
end

for i=1:numel(FFTAlignmentData)
    FFTAlignmentData(i).nopmean = nanmean(FFTAlignmentData(i).ordermat(:));
end

end

function [nfc] = directorFluctuationCorrelations(expt,FFTAlignmentData)
%% Format Data
% extract position information and convert to x & y matrices
[ax,ay] = xypos(FFTAlignmentData(1).pos);

% extract alignment angle and nematic order parameter data
ang = NaN([size(FFTAlignmentData(1).anglemat),size(FFTAlignmentData,2)]);
nop = NaN([size(FFTAlignmentData(1).ordermat),size(FFTAlignmentData,2)]); % nemtic order parameter
for i=1:size(FFTAlignmentData,2)
    ang(:,:,i) = FFTAlignmentData(i).anglemat;
    nop(:,:,i) = FFTAlignmentData(i).ordermat;
end

% Shift the angle range theta = [0,180] -> theta = [-90,90] & decompose
% into x & y vector components
adjang = 90;
pref = -1;
% Note - nop scales the alignment vectors by the nematic order parameter,
% which ranges from 0-1
[au,av] = alignmentvects(nop,pref.*ang,adjang);

%% Calculate Correlations
% Create plot directories 
corrtag = 'AlignmentFluctuations';
ctype = 'alignment';

% Iteratively calculate the correlations for each image frame
for i=1:size(au,3)
    disp(['Starting Frame ', num2str(i),'...']);
    [nfc(i)] = perpendicularFluctuationCorrelations(ax,ay,au(:,:,i),av(:,:,i),ax,ay,au(:,:,i),av(:,:,i),expt,i,corrtag);
    disp(['Finished Frame ', num2str(i),'!']);
end

end

function [dfc] = densityFluctuationCorrelations(expt,alignspec,stack,FFTAlignmentData)
%% calculate the density fluctuations
[dfluc] = calcDensityFluctuations(stack,expt,alignspec);

%% Load Density Fluctuation specific parameters if they exist, otherwise just use the same settings as the alignment analysis
parentPath = expt.path;
if exist(fullfile(parentPath,'DensityFlucParams.mat'))==2
    load(fullfile(parentPath,'DensityFlucParams.mat'));
else
    dflucspec = alignspec;
end

%% Format Alignment Vector Field Data
% extract position information and convert to x & y matrices
[ax,ay] = xypos(FFTAlignmentData(1).pos);

ang = NaN([size(FFTAlignmentData(1).anglemat),size(FFTAlignmentData,2)]);
nop = NaN([size(FFTAlignmentData(1).ordermat),size(FFTAlignmentData,2)]);
for i=1:size(FFTAlignmentData,2)
    ang(:,:,i) = FFTAlignmentData(i).anglemat;
    nop(:,:,i) = FFTAlignmentData(i).ordermat;
end
adjang = 90;
pref = -1;
[au,av] = alignmentvects(nop,pref.*ang,adjang);

%% Format Density Fluctuation Data
dx = dfluc(1).x;
dy = dfluc(1).y;
for i = 1:numel(dfluc)
    davg(:,:,i) = dfluc(i).mean;
end
du = davg;
dv = zeros(size(du)).*du;

%% Calculate Correlations
% Create plot directories 
corrtag = 'DensityFluctuations';
ctype = 'density';

for i=1:size(au,3)
    disp(['Starting Frame ', num2str(i),'...']);
    [dfc(i)] = perpendicularFluctuationCorrelations(dx,dy,du(:,:,i),dv(:,:,i),ax,ay,au(:,:,i),av(:,:,i),expt,i,corrtag);
    disp(['Finished Frame ', num2str(i),'!']);
end

end

function [fc] = perpendicularFluctuationCorrelations(x,y,u,v,ax,ay,au,av,expt,tidx,corrtag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x     - MxN position matrix in pixels
% y     - MxN position matrix in pixels
% u     - MxN matrix vector x-component - used to calculate correlations
% v     - MxN matrix vector y-component - used to calculate correlations
% ax    - x position matrix for alignment vectors
% ay    - y position matrix for alignment vectors
% au    - x-component alignment vector 
% av    - y-component alignment vector
% minr  - minimum distance for analysis
% maxr  - maximum distance for analysis
% exp   - structure containing experimental and analysis parameters
%
% Note: The alignment vectors are used to define the transverse and
% longitudinal axes for each grid position
%
% Computes the 1D correlation along a line thats transverse to direction of
% the orientation vector (nx,ny)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% check input params
if exist('perpcorr') == 0 || isempty(perpcorr), perpcorr = 0; end
if exist('interpdata') == 0 || isempty(interpdata), interpdata = 1; end

%% define the vector field for correlation
%Ensure your x & y positions are integer pixel values
xpix = round(ax);
ypix = round(ay);
%Get dimensions of vectors to be correlated
nrow = size(ax,1);
ncol = size(ax,2);

% analysis parameters
minrlen = 5;
rspace = 1;

%% Interpolate vector fields to the full image region
% generate x & y pixel location matrices
xs = 1:expt.stackdim(2);
ys = 1:expt.stackdim(1);
[xfull,yfull] = meshgrid(xs,ys);
intrpmth = 'cubic';

%Remove NaN's from vector fields for cubic and spline
%interpolation methods to prevent data truncation. If using a different
%interpolation method, the data will be preserved
%  Remove NaN's from vector field to be correlated
if interpdata == 1
    [ufull,vfull] = interpvfield(x,y,u,v,xfull,yfull,intrpmth);
else
    ufull = u;
    vfull = v;
end

%  Remove NaN's from alignment vector field
[aufull,avfull] = interpvfield(ax,ay,au,av,xfull,yfull,intrpmth);

%% Extract vector components from interpolated field at each grid position 
% Note:au_grid = au, av_grid = av for alignment vector correlations
for i=1:nrow % rows
    for j=1:ncol % cols
        ugrid(i,j) = ufull(ypix(i,j),xpix(i,j));
        vgrid(i,j) = vfull(ypix(i,j),xpix(i,j));
        augrid(i,j) = aufull(ypix(i,j),xpix(i,j));
        avgrid(i,j) = avfull(ypix(i,j),xpix(i,j));
    end
end

%Use registered vector components to generate a matrix of alignment vector
%angles at each grid point location for the vector field being correlation   
angle_grid = atan2(avgrid,augrid)*(180/pi);

%% Initialize structure for correlation data
% Fourier Correlation Structure
fc.c = repmat({NaN},[nrow,ncol]);
fc.clen = nan([nrow,ncol]);
fc.ang.longax = nan([nrow,ncol]);
fc.ang.transax = nan([nrow,ncol]);
fc0.c = repmat({NaN},[nrow,ncol]);

%% iterate over all vectors and perform correlation
I = nan(expt.stackdim(1),expt.stackdim(2));
n=min([2*expt.stackdim(1),2*expt.stackdim(2)]);
cont = 0;

for i=1:nrow
    tic
    for j=1:ncol
        if cont == 1
            if j ~= 1
                ii = i; jj = j-1;
            else
                ii = i-1; jj = ncol;
            end
            cont = 0;
        end
        disp([corrtag, ' -- Frame: ',num2str(tidx),' -- i(row) = ',num2str(i),' -- j(col) = ',num2str(j)])

        currang = angle_grid(i,j);
        % Ensure the angle measures stay on the same [-90,90] interval
        if (currang+90) <= 90
            perpang = currang + 90;
        else
            perpang = currang - 90;
        end
        %Get current position
        currpos = [xpix(i,j),ypix(i,j)];
        %Calculate the slope 'm' along the axis of the local alignment
        %vector as well as the slope 'mperp' of the perpendicular line
        m = tand(currang);
        mperp = tand(perpang);
        
        %Add the angles to the correlation data structures
        fc.ang.longax(i,j) = currang;
        fc.ang.transax(i,j) = perpang;
        
        %Find y-coordinates at which the perpendicular line intersects the
        %left-right edges of the image (at x=0, x=xmax)
        xend1 = [1,size(aufull,2)];
        yend1 = mperp.*(xend1-xpix(i,j))+ypix(i,j);
        
        %Find x-coordinates at which the perpendicular line intersects the
        %top-bottom edges of the image (at y=0, y=xmax)
        yend2 = [1,size(aufull,1)];
        xend2 = ((yend2-ypix(i,j))./mperp)+xpix(i,j);
        
        %Calculate which of the two lines (they are the same lines
        %extrapolated out to different boundaries) is shortest - use this
        %line as your transverse axis for correlation
        diff1 = abs(diff(yend1));
        diff2 = abs(diff(xend2));
        if isnan(diff1), diff1 = 0; end
        if isnan(diff2), diff2 = 0; end

        if diff1 == 0 && diff2 == 0
            % Skip remaining code if the current angle or position is NaN
            cont = 1;
            continue
        elseif diff1 < diff2
            xend = xend1;
            yend = yend1;
        else
            xend = xend2;
            yend = yend2;
        end
        
        %Determine indices of pixels along line
        [cx,cy,~] = improfile(I,xend,yend,n);
        cidx = [cx,cy];
        %Remove redundant pixel index pairs
        cidx = unique(round(cidx),'rows');
        %Remove any pixels outside of the image window - there shouldn't be
        %any if code above works properly
        ob = sum(cidx<=0,2);
        xob = cidx(:,1)>size(aufull,2);
        yob = cidx(:,2)>size(aufull,1);
        allob = ob + xob + yob;
        ob = allob>0;
        ib = ~ob;
        cidx = cidx(ib,:);
        
        %Preallocate vectors for x & y pixel coordinates along transverse
        %axis
        cposx = nan(size(cidx,1),1);
        cposy = nan(size(cidx,1),1);
        %Preallocate vectors for x & y components of vectors to be
        %correlated
        uperp = nan(size(cidx,1),1);
        vperp = nan(size(cidx,1),1);
        %For every position along the transverse axis, use the pixel index
        %to collect the corresponding x & y pixel coordinates as well as
        %the x & y vector components to be correlated at each position from
        %the full interpolated fields (where every pixel has a value)
        for k=1:size(cidx,1)
            cposx(k,1) = xfull(cidx(k,2),cidx(k,1));
            cposy(k,1) = yfull(cidx(k,2),cidx(k,1));
            uperp(k,1) = ufull(cidx(k,2),cidx(k,1));
            vperp(k,1) = vfull(cidx(k,2),cidx(k,1));
        end
        
        %Subtract the current grid position from the x & y pixel
        %coordinates along the axis to generate a vector of positions
        %relative to the current grid location
        %   x-y axes for new coordinates is congruent with the image axes
        perpx = cposx - currpos(1,1);
        perpy = cposy - currpos(1,2);
        
        %Calculate the relative distance from the current grid location to
        %each position along the transverse axis
        perpr = sqrt(perpx.^2+perpy.^2);

        %Calculate the relative angle between the current grid location and
        %each pixel along the transverse axis
        rang = atan2(perpy,perpx).*(180/pi);
        rperpang = rang-perpang;
        %Calculate the distance between the current grid location (zero
        %point) and each pixel along the transverse axis only (defined as
        %the current alignment vector direction + 90 degrees)
        rdist = perpr.*cosd(rperpang); 

        %Remove any nans from the x & y vector components along the
        %transverse axis - also remove nans from the positions along the
        %transverse axis (rdist)
        uu_nan = ~isnan(uperp);
        vv_nan = ~isnan(vperp);
        if ~isequal(uu_nan,vv_nan)
            error('nan issue!')
        end
        rdist = rdist(uu_nan);
        uperp = uperp(uu_nan);
        vperp = vperp(vv_nan); 
        
        if perpcorr == 1
            [uperp,vperp] = rotateCoordinateSpace(uperp,vperp,currang);
        end

        %Subtract the mean nematic director along the transverse axis from
        %each vector component - these are the fluctuations about the mean
        %nematic director
        u0perp = uperp-nanmean(uperp);
        v0perp = vperp-nanmean(vperp);   

        
        %Check to make sure our vectors weren't all nans
        if all(isnan(uperp)) == 1 || isempty(uperp) == 1 || size(uperp,1)==1
            cont = 1;
            continue
        end

        %Define equally spaced transverse axis positions for 1D
        %interpolation - spacing magnitude defined by rspace
        rInterp = [ceil(min(rdist)):rspace:floor(max(rdist))];
        rlen = numel(rInterp);

        %Check to ensure that our equally spaced transverse axis is longer
        %than our minimum length requirement (rlen)
        if rlen<minrlen
            cont = 1;
            continue
        end
        
        %Add the expected correlation lengths for the current grid position
        %to the correlation data structures - we will check to ensure that
        %the expected lengths match the actual lengths at the end
        if mod(rlen,2)>0
            fc.clen(i,j) = (rlen-1)/2;
        else
            fc.clen(i,j) = (rlen/2)-1;
        end
        
        %Interpolate equally spaced x & y fluctuation vector components
        %along the new transverse axis
        uInterp = interp1(rdist,u0perp,rInterp,'linear');
        vInterp = interp1(rdist,v0perp,rInterp,'linear');
        
        [radf] = fourierCorrelation1D(rInterp,uInterp,vInterp,expt);
        %Add outputs to correlation structures
        fc0.c{i,j} = radf.corr;
    end
    disp(['Finished Row ', num2str(i),' of ',num2str(size(x,1))]);
    toc
end

%% Check the correlation data for validity
[fc] = addnancells(fc0,[nrow,ncol]);

end

function [uff,vff] = interpvfield(x,y,u,v,xfull,yfull,intrpmth)
%Remove NaN's from alignment vector field for cubic and spline
%interpolation methods
if strcmp(intrpmth,'cubic') || strcmp(intrpmth,'spline')
    %Remove NaN's from vector field
    xnn = x; ynn = y; unn = u; vnn = v;
    [unn,urowidx,ucolidx] = delnans(unn);
    [vnn,vrowidx,vcolidx] = delnans(vnn);
    if ~isequal(urowidx,vrowidx), error('mismatch'); end
    if ~isequal(ucolidx,vcolidx), error('mismatch'); end
    xnn = xnn(urowidx,:); xnn = xnn(:,ucolidx);
    ynn = ynn(urowidx,:); ynn = ynn(:,ucolidx);
    %Interpolate without NaN's to prevent data truncation
    uff = interp2(xnn,ynn,unn,xfull,yfull,intrpmth);
    vff = interp2(xnn,ynn,vnn,xfull,yfull,intrpmth);
else
    %If using a different interpolation method, the data will be preserved
    uff = interp2(x,y,u,xfull,yfull,intrpmth);
    vff = interp2(x,y,v,xfull,yfull,intrpmth);
end
end

function [S] = addnancells(S,dim)
    nancells = cell(dim);
    nancells(:) = {nan};
    nanmat = nan(dim);
    names = fieldnames(S);
    for i=1:numel(names)
        tempcells = nancells;
        cfield = names{i};
        cval = S.(cfield);
        if iscell(cval)
            cval(cellfun(@isempty,cval(:))) = {nan};
            cdim = size(cval);
            tempcells(1:cdim(1),1:cdim(2)) = cval;
            S.(cfield) = tempcells;
        elseif isnumeric(cval)
            tempmat = nanmat;
            cdim = size(cval);
            tempmat(1:cdim(1),1:cdim(2)) = cval;
            S.(cfield) = tempmat;
        end
        
    end
end

function [S] = addemptycells(S,dim)
    nancells = cell(dim);
    names = fieldnames(S);
    for i=1:numel(names)
        temp = nancells;
        cfield = names{i};
        cval = S.(cfield);
        cdim = size(cval);
        temp(1:cdim(1),1:cdim(2)) = cval;
        S.(cfield) = temp;
    end
end

function [uperp,vperp] = rotateCoordinateSpace(uin,vin,currang)
%rotate the coordinate space to be perpendicualr and parallel to
%the current director
Rcc = [cosd(currang), -1.*sind(currang); sind(currang), cosd(currang)]; %counterclockwise
Rc = [cosd(currang), sind(currang); -1.*sind(currang), cosd(currang)]; % clockwise;
uv = [reshape(uin,[1 numel(uin)]);reshape(vin,[1 numel(vin)])];
for i=1:numel(uin)
    uvc(:,i) = Rc*uv(:,i);
end
para = uin'.*cosd(currang)-vin'.*sind(currang);
perp = uin'.*sind(currang)+vin'.*cosd(currang);
pp = [para;perp];
pperp = [zeros(size(para));perp];
for i=1:numel(uin)
    uvpcc(:,i) = Rcc*pperp(:,i);
end
uperp = uvpcc(1,:)';
vperp = uvpcc(2,:)';
end

function [radf] = fourierCorrelation1D(r,u,v,exp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates the spatial correlations of equally spaced
% vectors along a line/axis (1D) in real and fourier space
%                      
% All rights and permissions belong to:
% Ian Linsmeier
% ian.linsmeier@yale.edu
% Lab of Living Matter (Murrell Lab)
% 5/6/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check Inputs
if ~isequal(size(r),size(u)) || ~isequal(size(u),size(v))
   error('Dimensions of the inputs dont match!') 
end
if length(unique(diff(r))) ~= 1
    error('Vectors are not equally spaced')
end

%% Fourier Space Correlations
[radf] = fcorr(r,u,v);

end

function [acorr] = acorrvect(u,v)
%Calculates the unbiased autocorreltion of a 2D (x & y comp) vector series

[crx,lagx] = xcorr(u,'unbiased');
[cry,lagy] = xcorr(v,'unbiased');
if ~isequal(lagx,lagy), error('lagy code'); end

zidx = lagx>=0;
lagx = lagx(zidx);
lagy = lagy(zidx);
crx = crx(zidx);
cry = cry(zidx);

cr = crx + cry;
crn = cr./cr(1);

acorr.rpx = lagx;
acorr.corr = cr;
acorr.std = zeros(size(cr));
end

function [radf] = fcorr(r,u,v)
% check the size of each input
if ~all(isvector(r) && isvector(u) && isvector(v))
    error('inputs must be 1D vectors')
end
if ~isequal(size(r(:)),size(u(:))) || ~isequal(size(r(:)),size(v(:)))
    error('vectors must be the same length')
end
uFourier = fft(u);
vFourier = fft(v);
corrLen = numel(r);
deltaR = abs(r(2)-r(1));
deltaWaveNumber = 2*pi./deltaR;            
waveNumber = (mod(1/2+(0:(corrLen-1))/corrLen,1)-1/2).*deltaWaveNumber; 
fourierCorrelation = vFourier.*conj(vFourier) + uFourier.*conj(uFourier);
radf.corr = fourierCorrelation(waveNumber > 0);
end

function [dfluc] = calcDensityFluctuations(stack,expt,alignspec)
%% Define Threshold (range: [0,1]) below which image pixel values are set to NaN
thresh = 0.05;

%% Normalize stack and apply threshold
for i=1:size(stack,3)
    imadj(:,:,i) = imadjust(stack(:,:,i));
end
im = double(imadj);
imnorm = im./max(im(:));
imthresh = imnorm;
imthresh(find(imthresh<thresh)) = nan;
pctnan_total = nansum(isnan(imthresh(:)))./numel(imthresh(:))*100;
for i=1:size(stack,3)
    temp = isnan(imthresh(:,:,i));
    pctnan_im(i) = nansum(temp(:))./numel(temp(:))*100;
end

%% Calculate density fluctuations
if exist([pwd,'\DensityFlucParams.mat'])==2
    load([pwd,'\DensityFlucParams.mat']);
else
    dflucspec = alignspec;
end

[dfluc] = densityfluc(imthresh,dflucspec,expt);
end

function [dfluc] = densityfluc(im,dflucspec,expt)
winsz = dflucspec.winsize;
overlap = dflucspec.overlap;

winspace = winsz-round(winsz*overlap);
winrad = floor(winsz/2);

grid_row = winrad+1:winspace:expt.stackdim(1)-winrad;
grid_col = winrad+1:winspace:expt.stackdim(2)-winrad;
[xg,yg] = meshgrid(grid_col,grid_row);
numRows = length(grid_row);
numCols = length(grid_col);

strForm = sprintf('%%.%dd',length(num2str(numRows)));
fprintf('Starting row ')

for f=1:expt.stackdim(3)
    for i = 1:numRows
        procStr = sprintf(strForm,i);
        fprintf(1,[procStr,' of ',num2str(numRows),' - Frame: ',num2str(f),'\n']);
        for j = 1:numCols
            window = im(grid_row(i)-winrad:grid_row(i)+winrad,grid_col(j)-winrad:grid_col(j)+winrad,f);
            winmean(i,j,f) = nanmean(window(:));
            winstd(i,j,f) = nanstd(window(:));
            winsqrt(i,j,f) = sqrt(nanmean(window(:)));
        end
    end
    dfluc(f).winsz = winsz;
    dfluc(f).overlap = overlap;
    dfluc(f).winspace = winspace;
    dfluc(f).winrad = winrad;
    dfluc(f).x = xg;
    dfluc(f).y = yg;
    dfluc(f).mean = winmean(:,:,f);
    dfluc(f).std = winstd(:,:,f);
    dfluc(f).nstd = winstd(:,:,f)./winmean(:,:,f);
    dfluc(f).sqrt = winsqrt(:,:,f);
end
end

function sanitycheck1(x,y,u,v,ax,ay,au,av,xpix,ypix,ufull,vfull,aufull,avfull,ugrid,vgrid,augrid,avgrid,angle,angle_grid,angle_grid2)
    if isequaln(x,xpix) == 0, error('WTF!'); end
    if isequaln(y,ypix) == 0, error('WTF!'); end
    if isequaln(x,ax) == 0, error('WTF!'); end
    if isequaln(y,ay) == 0, error('WTF!'); end
    if isequaln(u,au) == 0, error('WTF!'); end
    if isequaln(v,av) == 0, error('WTF!'); end
    if isequaln(ufull,aufull) == 0, error('WTF!'); end
    if isequaln(vfull,avfull) == 0, error('WTF!'); end
    if isequaln(ugrid,augrid) == 0, error('WTF!'); end
    if isequaln(vgrid,avgrid) == 0, error('WTF!'); end
    if isequaln(angle_grid,angle_grid2) == 0, error('WTF!'); end
    % When we register the alignment vector field with itself we should get the
    % same matrix back when we extract the data at the grid locations (au =
    % au_grid, av = av_grid)
    if isequaln(au,augrid) == 0, error('WTF!'); end
    if isequaln(av,avgrid) == 0, error('WTF!'); end
    if isequaln(angle_grid,angle) == 0, error('WTF!'); end
end

function sanitycheck2(i,j,currang,perpang,m,mperp,xend,yend,expt)
    disp(['i(row) = ',num2str(i),'; j(col) = ',num2str(j),';'])
    %             disp(['row = ',num2str(i),'; col = ',num2str(j),';'])
    disp(['currang = ',num2str(currang)])
    disp(['perpang = ',num2str(perpang)])
    disp(['m = ',num2str(m)])
    disp(['mperp = ',num2str(mperp)])
    disp(['xend = [',num2str(xend(1)),', ',num2str(xend(2)),']'])
    disp(['yend = [',num2str(yend(1)),', ',num2str(yend(2)),']'])
    if isequal(xend,[1,expt.stackdim(2)])
        disp(sprintf(['Transverse axis intersects the image edges at:\nLeft Edge Intersection: (x = ',num2str(xend(1)),', y = ',num2str(yend(1)),')\nRight Edge Intersection: (x = ',num2str(xend(2)),', y = ',num2str(yend(2)),')']))
    else
        disp(sprintf(['Transverse axis intersects the image edges at:\nTop Edge Intersection: (x = ',num2str(xend(1)),', y = ',num2str(yend(1)),')\nBottom Edge Intersection: (x = ',num2str(xend(2)),', y = ',num2str(yend(2)),')']))
    end     
end

function sanitycheck3(fc,rc,ang)
%Check the lengths of the correlations in fourier space
fclen0 = cellfun(@numel,fc.c);
fcnan = cell2mat(cellfun(@(x) all(isnan(x(:))),fc.c,'uniform',0));
fclen = fclen0-fcnan;

fclen2 = fc.clen;
fclen2(isnan(fclen2)) = 0;

if ~isequal(fclen,fclen2)
    error('Something went terribly wrong...')
end

%Check the lengths of the correlations in real space
rclen0 = cellfun(@numel,rc.c);
rcnan = cell2mat(cellfun(@(x) all(isnan(x(:))),rc.c,'uniform',0));
rclen = rclen0-rcnan;

rclen2 = rc.clen;
rclen2(isnan(rclen2)) = 0;

if ~isequal(rclen,rclen2)
    error('Something went terribly wrong...')
end

if ~isequaln(fc.ang.longax,ang) || ~isequaln(rc.ang.longax,ang)
    error('Something went terribly wrong...')
end
end

function [stack,stackInfo] = stackread(fullpath)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program reads in a multiplaned TIFF file
% into a single variable
%
% 11/12/09 Patrick Oakes
% poakes@gmail.com
%
% Updated
% 5/26/16 - Ian Linsmeier
% - Changed function name
% - Compatible with 8,16, & 32-bit images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('fullpath') || isempty(fullpath)
    [file,path] = uigetfile({'*.tif';'*.png';'*.jpg';'*.pdf';'*.gif'});
    if isequal(path, 0) || isequal(file, 0), return; end
    fullpath = fullfile(path,file);
end

stackInfo = imfinfo(fullpath);
dim = size(stackInfo);
if dim(1) > dim(2)
    nFrames = dim(1);
else
    nFrames = dim(2);
end
stack = zeros(stackInfo(1).Height,stackInfo(1).Width,nFrames);

if double(stackInfo(1).BitDepth) == 24
    stack = uint8(zeros(stackInfo(1).Height,stackInfo(1).Width,3,nFrames));
end

if nFrames == 1
    stack = imread(fullpath); 
else
    for i = 1:nFrames
        if double(stackInfo(1).BitDepth) == 24
            stack(:,:,:,i) = imread(fullpath,i);
        else
            stack(:,:,i) = imread(fullpath,i);
        end
    end
end

if double(stackInfo(1).BitDepth) == 8
    stack = uint8(stack);
elseif double(stackInfo(1).BitDepth) == 16
    stack = uint16(stack);
elseif double(stackInfo(1).BitDepth) == 32
    stack = uint32(stack);
end
end

function [u,v,angshift] = alignmentvects(order,angle,adjang)
% angle - vector or matrix of angle measures in degrees
% order - vector  or matrix of nematic order parameter magnitude
% adjang - angle measure (degrees) by which to rotate the angle vector

% % Check to ensure inputs are column or row vectors and their dimensions
% % match
% if sum(ismember(size(angle),1))==0 || sum(ismember(size(order),1))==0 
%     error('Inputs must be Nx1 or 1xN vectors')
% end
% if size(angle)~=size(order)
%     error('Inputs must have the same dimensions')
% end
% 
% Check to see if angles are in radians - Range is less than [-2*pi,+2*pi]

if range(angle(:))<=4*pi
    disp('WARNING: angle values must be in degrees!')
    disp('Current angle range seems to low')
end

% angle = mod(wrapTo360(angle+adjang),180);
% angle = mod(angle+adjang,180);
angshift = angle+adjang;

u = order.*cosd(angshift);
v = order.*sind(angshift);

end

function [A,rowidx,colidx] = delnans(A)
% Removes rows and columns from a 2D matrix A that are all nans

dim = size(A);
if length(dim)>2
    error('The matrix must be 2D!')
end
rowidx = ~all(isnan(A),2);
colidx = ~all(isnan(A));
A = A(rowidx,:); % for nan - rows
A = A(:,colidx);   % for nan - columns
end

function [xx,yy] = xypos(u,v)
% Generates a matrix of x & y positions from a single x-y position vector
% u, or individual x & y position vectors u & v, respectively 

if exist('v')==0 || isempty(v)==1
    if size(u,1) < size(u,2)
        u = u';
    end
    if size(u,2) == 2
        v = u(:,2);
        u = u(:,1);
    else
        error('Unexpected Position Vector Size')
    end
end

x = unique(u);
y = unique(v);
[xx,yy] = meshgrid(x,y);
end




