function [mask, info] = actinFilamentMask(I, varargin)
%ACTINFILAMENTMASK  Segment F-actin filaments from a single phalloidin/TxRed image.
%
%   Targets LONG, COHERENT FILAMENTS — cortical fibers running parallel to the
%   cell boundary and stress fibers projecting toward it — which are most
%   prominent at the cell periphery. It is not trying to trace a mesh, and the
%   defaults actively suppress fine cytoplasmic texture.
%
%   MASK = ACTINFILAMENTMASK(I) returns a logical mask for one 2-D image I
%   (uint16 SUM z-projection, linear intensities).
%
%   [MASK, INFO] = ACTINFILAMENTMASK(...) also returns intermediates and QC.
%
%   Pipeline, and why each step is the way it is:
%     1. Mild denoise.
%     2. Provisional Otsu mask -> used ONLY to site the noise estimate. It
%        lands on bright perinuclear cytoplasm, which is a fine place to
%        measure noise and a terrible place to look for filaments.
%     3. Noise sigma from the RAW image via a Laplacian kernel, inside the
%        provisional mask. Zero response to linear ramps, so structure drops
%        out; restricted to cells because a background clipped at zero carries
%        no noise and would drag a whole-field estimate down.
%     4. PERMISSIVE cell mask: smoothed > CellThreshold x noiseSigma. Otsu is
%        NOT used here. Measured on B03_S03: Otsu kept 26% of the field while
%        cells covered 86%, and only 3.6% of the 0-2 µm peripheral zone fell
%        inside it — so the fibers of interest were outside the search region
%        entirely. Otsu separates bright cytoplasm from everything else; it
%        does not separate cells from background.
%     5. White top-hat -> removes background gradient and the diffuse actin
%        pool, keeps thin structures.
%     6. LOCAL CONTRAST NORMALIZATION -> the top-hat is converted to a local
%        z-score using statistics gathered only from cell pixels. Without this,
%        absolute contrast rises toward the bright interior (mean top-hat 121
%        at the periphery vs 270 in the interior on S03) and the detector
%        prefers perinuclear texture over real peripheral fibers. In z, the
%        ranking inverts and matches the biology: ridge-like signal is densest
%        at the periphery.
%     7. FIBERMETRIC on the z-image. StructureSensitivity is therefore in
%        DIMENSIONLESS z UNITS, not intensity units — it no longer depends on
%        how bright the field is.
%     8. Hysteresis on the ridge response with a z floor.
%     9. LENGTH FILTER -> keeps long objects, drops short ones. This is what
%        enforces "filaments, not mesh".
%
%   Name-value arguments:
%     'PixelSize'        µm/px. Default 0.35.
%     'ValidMask'        Analysable region (e.g. excluding the deconvolution
%                        border). Default [] (whole field).
%     'CellThreshold'    Cell mask cut, in multiples of noise sigma on the
%                        smoothed image ABOVE the estimated background
%                        baseline (see 'BackgroundPercentile'). Default 10.
%                        Lower to include fainter lamella, raise if
%                        background creeps in. Check INFO.CellFraction
%                        against what you see.
%     'BackgroundPercentile' Percentile of the smoothed, valid-region image
%                        used to estimate the background baseline that
%                        gets subtracted before the CellThreshold cut.
%                        Default 10. Needed because a comparison against
%                        raw S alone assumes background sits near zero --
%                        true only if upstream processing left it that way
%                        (e.g. restoreOffset=false). If your SUM_linear
%                        files carry a real, non-zero background baseline
%                        (restoreOffset=true), this must be estimated and
%                        subtracted or CellFraction will come out near 1.0.
%     'FilamentWidth'    [min max] µm. Default [0.7 2.5]. The lower bound
%                        matters: below ~0.7 µm the scales reach cytoplasmic
%                        granularity and the mask fills with texture.
%     'TopHatRadius'     µm. Default 4.
%     'NormWindow'       µm. Window for local contrast normalization.
%                        Default 10. Should exceed filament spacing but stay
%                        below cell size.
%     'StructSens'       StructureSensitivity for FIBERMETRIC, in z units.
%                        Default 2.
%     'LocalSdCap'       Local SD is capped at this multiple of its field-wide
%                        median before normalizing. Default 2. Without it,
%                        bright bundles inflate their own local SD and get
%                        divided back below threshold. Raise toward Inf to
%                        disable; lower to rescue bundles harder.
%     'RidgeThreshold'   High hysteresis threshold on the ridge response.
%                        Default 0.04. Low threshold is a quarter of it.
%     'ZFloor'           Local z a pixel must clear regardless of ridge
%                        response. Default 1.5.
%     'MinLength'        Major axis below this (µm) is dropped. Default 3 —
%                        this is the filaments-not-mesh knob.
%     'MinArea'          µm^2. Default 0.5.
%     'DenoiseSigma'     µm. Default 0.2. 0 disables.
%     'CellMask'         Supply your own to skip steps 2 and 4.
%
%   Tuning: LOOK AT THE OVERLAYS. Every summary statistic tried on this data
%   (mean width, object count, coverage) failed to separate a good mask from a
%   bad one. INFO's numbers are for provenance and for checking consistency
%   ACROSS fields once a setting is chosen by eye — not for choosing it.
%
%   Requires the Image Processing Toolbox.
%
%   Example:
%     I  = imread('B03_S03_TxRed_SUM_linear.tif');
%     vm = imread('B03_S01_TxRed_validmask.tif') > 0;
%     [m, info] = actinFilamentMask(I, 'ValidMask', vm);
%     imshow(imoverlay(imadjust(I), m, [0 1 0]))

%% ---- parse inputs -------------------------------------------------------
p = inputParser;
p.addRequired('I', @(x) isnumeric(x) && ismatrix(x));
p.addParameter('PixelSize', 0.35, @(x) isscalar(x) && x > 0);
p.addParameter('ValidMask', [], @(x) isempty(x) || islogical(x));
p.addParameter('CellThreshold', 10, @(x) isscalar(x) && x > 0);
p.addParameter('BackgroundPercentile', 10, @(x) isscalar(x) && x >= 0 && x <= 100);
p.addParameter('FilamentWidth', [0.7 2.5], @(x) numel(x) == 2 && all(x > 0));
p.addParameter('TopHatRadius', 4, @(x) isscalar(x) && x > 0);
p.addParameter('NormWindow', 10, @(x) isscalar(x) && x > 0);
p.addParameter('StructSens', 2, @(x) isscalar(x) && x > 0);
p.addParameter('LocalSdCap', 2, @(x) isscalar(x) && x > 0);
p.addParameter('RidgeThreshold', 0.04, @(x) isscalar(x) && x > 0 && x < 1);
p.addParameter('ZFloor', 1.5, @(x) isscalar(x) && x >= 0);
p.addParameter('MinLength', 3, @(x) isscalar(x) && x >= 0);
p.addParameter('MinArea', 0.5, @(x) isscalar(x) && x >= 0);
p.addParameter('DenoiseSigma', 0.2, @(x) isscalar(x) && x >= 0);
p.addParameter('CellMask', [], @(x) isempty(x) || islogical(x));
p.parse(I, varargin{:});
opt = p.Results;

px  = opt.PixelSize;
Raw = double(I);

if isempty(opt.ValidMask)
    validMask = true(size(Raw));
else
    validMask = opt.ValidMask;
    if ~isequal(size(validMask), size(Raw))
        error('actinFilamentMask:ValidMaskSize', 'ValidMask must match size(I).');
    end
end

%% ---- 1. denoise --------------------------------------------------------
sigPx = opt.DenoiseSigma / px;
if sigPx > 0, Id = imgaussfilt(Raw, sigPx); else, Id = Raw; end
S = imgaussfilt(Id, 3 / px);              % ~3 µm blur: cells, not fibers

%% ---- 2-3. noise sigma, sited inside bright cytoplasm -------------------
Sn   = mat2gray(S, prctile(S(validMask), [0.1 99.9]));
prov = (Sn > graythresh(Sn(validMask))) & validMask;   % provisional, for noise only
if ~any(prov(:))
    error('actinFilamentMask:NoForeground', 'No foreground found — blank image?');
end

K  = [1 -2 1; -2 4 -2; 1 -2 1];          % zero response to linear ramps; norm 6
L  = imfilter(Raw, K, 'symmetric');
Lv = L(prov);
noiseSigma = 1.4826 * median(abs(Lv - median(Lv))) / 6;
if noiseSigma <= 0, noiseSigma = std(Lv) / 6; end
if noiseSigma <= 0
    error('actinFilamentMask:ZeroNoise', 'Noise sigma estimated as zero.');
end

%% ---- 4. permissive cell mask -------------------------------------------
if isempty(opt.CellMask)
    % S > CellThreshold*noiseSigma alone silently assumes background sits
    % near zero. That was true while the upstream deconvolution pipeline
    % flattened background to (or near) zero, but is no longer true once
    % restoreOffset correctly restores the camera's real, non-zero
    % baseline -- comparing straight against S would then call the whole
    % field "cell". Subtract a robust background baseline first so the
    % comparison is relative to actual background, not to zero.
    backgroundLevel = prctile(S(validMask), opt.BackgroundPercentile);
    cellMask = ((S - backgroundLevel) > opt.CellThreshold * noiseSigma) & validMask;
    cellMask = imfill(cellMask, 'holes');
    cellMask = imopen(cellMask, strel('disk', round(1 / px)));
    cellMask = bwareaopen(cellMask, round(50 / px^2));   % drop <50 µm^2 debris
else
    cellMask = opt.CellMask & validMask;
    backgroundLevel = NaN; % not estimated when a custom CellMask is supplied
end
if ~any(cellMask(:))
    error('actinFilamentMask:EmptyCellMask', ...
          'Cell mask is empty — CellThreshold may be too high.');
end

%% ---- 5. top-hat --------------------------------------------------------
rTH = max(2, round(opt.TopHatRadius / px));
TH  = imtophat(Id, strel('disk', rTH));

%% ---- 6. local contrast normalization -----------------------------------
% Statistics are gathered from cell pixels only (normalized convolution).
% Letting the dark exterior into the window would depress the local mean near
% every cell edge and paint a false bright rim around all of them.
w  = 2 * floor(opt.NormWindow / px / 2) + 1;
h  = ones(w) / w^2;
C  = double(cellMask);
n1 = imfilter(C,           h, 'symmetric');
m1 = imfilter(TH .* C,     h, 'symmetric') ./ max(n1, eps);
m2 = imfilter(TH.^2 .* C,  h, 'symmetric') ./ max(n1, eps);
localSd = sqrt(max(m2 - m1.^2, 0));

% Cap the local SD. A bright bundle IS its own neighbourhood, so it inflates
% the very statistic used to normalize it and suppresses itself. Measured at
% x=911,y=637 in B03_S03: top-hat 3294 (109x noise), localSd 948 (99.7th pct
% of the field), Z = 1.24 -> silently dropped. Capping at 2x the field median
% lifts that bundle to Z = 3.13 while moving the fraction of cell pixels above
% ZFloor by only 7.4% -> 7.9%. It bounds how far any structure can suppress
% itself without touching the dim peripheral fibers.
medSd = median(localSd(cellMask));
sdCap = opt.LocalSdCap * medSd;
denom = max(min(localSd, sdCap), noiseSigma);

Z = (TH - m1) ./ denom;
Z(~cellMask) = 0;

%% ---- 7. ridge response on the z-image ----------------------------------
wPx       = sort(opt.FilamentWidth) / px;
thickness = unique(max(1, round(wPx(1)) : round(wPx(2))));

V  = fibermetric(Z, thickness, 'ObjectPolarity', 'bright', ...
                 'StructureSensitivity', opt.StructSens);
Vn = min(max(V, 0), 1);

%% ---- 8. hysteresis -----------------------------------------------------
tHi = opt.RidgeThreshold;
tLo = 0.25 * tHi;

seed = (Vn > tHi) & (Z > opt.ZFloor)       & cellMask;
grow = (Vn > tLo) & (Z > 0.5 * opt.ZFloor) & cellMask;
mask = imreconstruct(seed, grow);

%% ---- 9. size and length filtering --------------------------------------
mask = mask & validMask;
mask = bwareaopen(mask, max(1, round(opt.MinArea / px^2)));

if opt.MinLength > 0 && any(mask(:))
    mask = bwpropfilt(mask, 'MajorAxisLength', [opt.MinLength / px, Inf]);
end
mask = bwmorph(mask, 'clean');

%% ---- info --------------------------------------------------------------
if nargout > 1
    skel = bwskel(mask);
    info = struct( ...
        'ValidMask',        validMask, ...
        'CellMask',         cellMask, ...
        'ProvisionalMask',  prov, ...
        'TopHat',           TH, ...
        'Zscore',           Z, ...
        'MedianLocalSd',    medSd, ...
        'SdCap',            sdCap, ...
        'Ridge',            Vn, ...
        'ThicknessPx',      thickness, ...
        'NoiseSigma',       noiseSigma, ...
        'BackgroundLevel',  backgroundLevel, ...
        'CellFraction',     nnz(cellMask) / nnz(validMask), ...
        'ProvFraction',     nnz(prov)     / nnz(validMask), ...
        'SkeletonLength_um', nnz(skel) * px, ...
        'Coverage',         nnz(mask) / nnz(validMask), ...
        'CoverageInCells',  nnz(mask) / nnz(cellMask));
end
end
