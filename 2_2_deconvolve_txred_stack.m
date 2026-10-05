function R = deconvolve_txred_stack(fileList, opts)
%DECONVOLVE_TXRED_STACK  Deconvolve one 9-slice TxRed (F-actin) z-stack.
%
%   R = deconvolve_txred_stack(fileList, opts)
%
% INPUT
%   fileList  cellstr/string array of the slice paths for ONE field.
%             Any order - sorted internally by the Z token in the filename.
%   opts      struct:
%     REQUIRED
%       .NA            objective numerical aperture        (e.g. 0.80)
%       .n_medium      immersion refractive index          (e.g. 1.00 air)
%       .pixelSize     micron/pixel                        (e.g. 0.35)
%       .zStep         micron between slices               (e.g. 2.00)
%       .cameraOffset  ONE fixed ADU value for the TxRed channel. See NOTE 2.
%     OPTIONAL (defaults)
%       .scaleFactor      1/9    constant; see NOTE 3
%       .numIter          15     RL iterations
%       .doDeconvolution  true   false -> raw sum projection (still linear)
%       .psfNZ            5      PSF depth in planes (odd)
%       .psfOversample    5      odd -> PSF peak centred in a pixel
%       .edfSmooth        3      px smoothing of the EDF focus map
%       .restoreOffset    false  true -> background back near raw ADU (AFFINE)
%       .outBase          ''     non-empty -> also writes <outBase>_*.tif/png
%
% OUTPUT
%   R.sum   uint16  LINEAR, single-slice-like scale -> QUANTIFY FROM THIS
%   R.seg   uint16  EDF all-in-focus, NON-linear -> tracing/segmentation ONLY
%   R.mask  logical usable interior (border ringing excluded)
%   R.info  struct  one row of diagnostics for your own manifest
%
% EXAMPLE LOOP (yours to adapt - discovers fields, no hardcoded lists)
%   opts = struct('NA',0.80,'n_medium',1.00,'pixelSize',0.35,'zStep',2.00, ...
%                 'cameraOffset',601);      % <- from a DARK FRAME, see NOTE 2
%   inFolder  = pwd;
%   outFolder = fullfile(inFolder,'processed');
%   if ~exist(outFolder,'dir'), mkdir(outFolder); end
%
%   d   = dir(fullfile(inFolder,'Well*_ChannelTxRed_Seq*_Z*.tif'));
%   tok = regexp({d.name}, ...
%       '^Well([^_]+)_Channel([^_]+)_Seq(\d+)_S(\d+)_Z(\d+)\.tif$','tokens','once');
%   keep = ~cellfun(@isempty, tok);   d = d(keep);   tok = tok(keep);
%   key  = cellfun(@(t) sprintf('%s_S%s', t{1}, t{4}), tok, 'UniformOutput', false);
%
%   uk = unique(key);   rows = [];
%   for i = 1:numel(uk)
%       sel = strcmp(key, uk{i});
%       f   = fullfile({d(sel).folder}, {d(sel).name});
%       if numel(f) ~= 9      % also catches a field acquired twice (two Seq)
%           warning('%s: %d slices, expected 9 - skipped.', uk{i}, numel(f));
%           continue
%       end
%       opts.outBase = fullfile(outFolder, [uk{i} '_TxRed']);
%       R = deconvolve_txred_stack(f, opts);
%       rows = [rows; R.info];
%   end
%   writetable(struct2table(rows), fullfile(outFolder,'txred_manifest.csv'));
%
%   % then confirm consistency across every stack you intend to compare:
%   %   all offsetUsed / numIter / scaleFactor identical, fluxRatio ~ 1.000,
%   %   nSaturatedRaw and nClippedOut all zero.
%
% -------------------------------------------------------------------------
% NOTES - these decide whether your numbers mean anything
% -------------------------------------------------------------------------
% 1. AT NA 0.8 / 0.35 um px / 2 um z-step YOUR SYSTEM IS UNDERSAMPLED.
%    TxRed: 1.12 px per lateral FWHM, 0.85 slices per axial FWHM (Nyquist >=2).
%    Deconvolution CANNOT recover lateral resolution here and will NOT make
%    filaments thinner - the image is pixel-limited, not diffraction-limited.
%    What you gain is out-of-focus haze removal, i.e. contrast.
%    A Gaussian PSF is INVALID at this sampling (sigma_z ~ 0.3 slices = a delta
%    function; it models no haze and measured WORSE than raw on your data).
%    This function builds a scalar-diffraction widefield PSF from the pupil
%    function, which reproduces the defocus cone (~15 px radius at 4 um).
%
% 2. cameraOffset MUST BE ONE FIXED NUMBER, AND FOR TxRed THIS MATTERS MOST.
%    A per-image estimate is NOT safe on this channel. Measured on your data
%    the per-image estimate is 601 (B03/S01) vs 634 (C03/S02) - 33 ADU, ~4% of
%    signal - because an actin image has NO empty background population: the
%    histogram mode holds only ~0.1-0.2% of pixels and just tracks diffuse
%    haze, i.e. the SAMPLE, not the camera. Let it float per well and a
%    treatment that changes diffuse actin gets silently subtracted away.
%    Also note DAPI floors near 414 while TxRed floors near 601 on the same
%    camera: a shared black level cannot differ by ~190 ADU, so most of the
%    TxRed floor is REAL diffuse background (media autofluorescence, unbound
%    phalloidin, haze), not an offset. Subtracting it is a background-
%    subtraction decision - make it deliberately and use one value everywhere.
%    Correct source: a dark frame (no illumination, same exposure/gain).
%
% 3. scaleFactor MUST BE CONSTANT ACROSS EVERY STACK YOU COMPARE.
%    Default 1/9 turns the z-sum into a z-average so values land in
%    single-slice range. Deliberately a fixed constant, NOT 1/nZ: a stack with
%    fewer slices genuinely holds less signal and 1/nZ would hide that.
%    Recover absolute summed flux with flux = double(R.sum)/R.info.scaleFactor.
%    Deconvolution concentrates light into filament cores, so watch for the
%    clipping warning - if it fires, lower scaleFactor AND rerun every stack.
%
% 4. FOR PURE PHOTOMETRY, DECONVOLUTION IS OPTIONAL. Integrated intensity is
%    conserved either way (measured flux ratio 1.0000), so a sum projection of
%    the offset-subtracted raw stack is already correct for total actin. Run
%    once with doDeconvolution=false and compare. If your numbers move,
%    deconvolution is redistributing signal. Its real value is SEGMENTATION.
%
% 5. R.sum is LINEAR; R.seg is NOT. R.seg is an EDF (extended depth of field)
%    blend: each pixel is a sharpness-weighted mean across z. The weights vary
%    spatially, so R.seg values are NOT proportional to fluorophore. It is the
%    better image for TRACING filaments and the wrong one for measuring them.
%
% 6. FILAMENT WIDTH IS NOT A MEASUREMENT HERE. Richardson-Lucy narrows
%    structures by an iteration-dependent amount and never converges, so width
%    and peak intensity drift with numIter. Never compare widths across
%    different numIter. Over-iterating also splits one thick bundle into
%    several thin ones, inflating filament counts convincingly. Keep numIter
%    identical everywhere and sanity-check against the raw best-focus slice.
%
% Requires: Image Processing Toolbox.
% -------------------------------------------------------------------------

CHANNEL   = 'TxRed';
LAMBDA_EM = 615;     % Texas Red emission, nm - longer than DAPI => wider PSF
SEG_LABEL = 'EDF';   % thin filaments: sharpness-weighted focus stacking is
                     % ideal (unlike nuclei, where it speckles smooth interiors)

% PSF cache: the PSF depends only on the optics, so in your loop it is built
% once and reused for every stack thereafter. Must be declared at function
% top-level, not inside a conditional.
persistent cachedPSF cachedKey

%% ---- 0. Options ------------------------------------------------------
if nargin < 2, opts = struct(); end
req = {'NA','n_medium','pixelSize','zStep','cameraOffset'};
for i = 1:numel(req)
    if ~isfield(opts, req{i}) || isempty(opts.(req{i}))
        error('deconvolve_txred_stack:missingOpt', ...
            ['opts.%s is required.\n' ...
             'cameraOffset in particular has no ''auto'' default by design: on ' ...
             'TxRed a per-image estimate tracks diffuse haze rather than the ' ...
             'camera (601 vs 634 across your two wells), so it would normalise ' ...
             'each condition separately and invalidate the comparison.'], req{i});
    end
end
opts = setdef(opts, 'scaleFactor',     1/9);
opts = setdef(opts, 'numIter',         15);
opts = setdef(opts, 'doDeconvolution', true);
opts = setdef(opts, 'psfNZ',           5);
opts = setdef(opts, 'psfOversample',   5);
opts = setdef(opts, 'edfSmooth',       3);
opts = setdef(opts, 'restoreOffset',   false);
opts = setdef(opts, 'outBase',         '');

%% ---- 1. Order and load the stack (raw ADU, no scaling) ---------------
fileList = cellstr(fileList);
tok = regexp(fileList, '[Zz](\d+)', 'tokens', 'once');
if all(~cellfun(@isempty, tok))
    [~, ord] = sort(cellfun(@(c) str2double(c{1}), tok));
    fileList = fileList(ord);
else
    warning('deconvolve_txred_stack:noZToken', ...
            'No Z token in every filename - using the order given.');
end
nZ = numel(fileList);

info1 = imfinfo(fileList{1});
raw = zeros(info1(1).Height, info1(1).Width, nZ, 'double');
for k = 1:nZ
    a = double(imread(fileList{k}));
    if ~isequal(size(a), [info1(1).Height info1(1).Width])
        error('Slice %d (%s) size differs from slice 1.', k, fileList{k});
    end
    raw(:,:,k) = a;
end
rawMin = min(raw(:));  rawMax = max(raw(:));
nSat = sum(raw(:) >= 65535);
if nSat > 0
    warning('deconvolve_txred_stack:saturated', ...
        ['%d saturated px at 65535 - integrated intensity is UNDERESTIMATED ' ...
         'for those filaments.'], nSat);
end

%% ---- 2. Sampling diagnostics ----------------------------------------
lam_um    = LAMBDA_EM/1000;
FWHMxy    = 0.51 * lam_um / opts.NA;
FWHMz     = 1.77 * opts.n_medium * lam_um / opts.NA^2;
pxPerFWHM = FWHMxy / opts.pixelSize;
szPerFWHM = FWHMz  / opts.zStep;
if pxPerFWHM < 2
    warning('deconvolve_txred_stack:undersampledXY', ...
        ['Laterally undersampled %.1fx (%.2f px per FWHM). No lateral ' ...
         'resolution recovery is possible - expect haze removal only, not ' ...
         'thinner filaments.'], opts.pixelSize/(FWHMxy/2), pxPerFWHM);
end

%% ---- 3. Background: fixed constant subtraction (AFFINE = linear-safe)
offsetUsed = double(opts.cameraOffset);
imgDiff = raw - offsetUsed;
nClamped = sum(imgDiff(:) < 0);
% RL needs non-negative input, but max(x,0) alone throws away the actual
% raw-vs-offset fluctuation for every pixel that happened to read below
% the offset (which real read noise does ~half the time near a dark
% level). Keep that discarded part: since max(x,0)+min(x,0) == x exactly,
% adding negResidual back after deconvolution (see step 6) reconstructs
% the true raw value in background regions instead of flattening them to
% one constant.
negResidual = min(imgDiff, 0);
img = max(imgDiff, 0);   % RL needs non-negativity; only sub-floor background hit

%% ---- 4. PSF (cached: built once per unique optical parameter set) ----
psf = [];
if opts.doDeconvolution
    psfNZ = opts.psfNZ;
    if mod(psfNZ,2) == 0, psfNZ = psfNZ - 1; end
    psfNZ = min(psfNZ, nZ);

    persistent_key = [opts.NA opts.n_medium opts.pixelSize opts.zStep psfNZ ...
           opts.psfOversample LAMBDA_EM];
    if isempty(cachedKey) || ~isequal(persistent_key, cachedKey)
        cachedPSF = widefieldPSF(LAMBDA_EM, opts.NA, opts.n_medium, ...
                                 opts.pixelSize, opts.zStep, psfNZ, opts.psfOversample);
        cachedKey = persistent_key;
    end
    psf = cachedPSF;    % identical PSF for every stack in your loop
else
    psfNZ = 0;
end

%% ---- 5. Deconvolve ---------------------------------------------------
if opts.doDeconvolution
    rXY = (size(psf,1)-1)/2;  rZ = (size(psf,3)-1)/2;
    padXY = 2*rXY;  padZ = min(rZ, nZ-1);
    imgPad   = padarray(img, [padXY padXY padZ], 'symmetric', 'both');
    deconPad = deconvlucy(imgPad, psf, opts.numIter, 0);   % dampar=0 -> unbiased
    decon = deconPad(padXY+1:end-padXY, padXY+1:end-padXY, padZ+1:end-padZ);
    decon = max(decon, 0);
    clear imgPad deconPad
    fluxRatio = sum(decon(:)) / max(sum(img(:)), eps);
    if abs(fluxRatio - 1) > 0.05
        warning('deconvolve_txred_stack:flux', ...
            'Flux ratio %.4f drifted >5%% from 1 - investigate before quantifying.', ...
            fluxRatio);
    end
else
    decon = img;  rXY = 2;  fluxRatio = 1;
end

%% ---- 6a. LINEAR photometric image (QUANTIFY FROM THIS) ---------------
% Integrated intensity over z is proportional to total actin per (x,y). Sum is
% linear; scaleFactor is a constant multiply, so proportionality survives.
if opts.restoreOffset
    % Add back the exact clamped residual (not just the flat offset) so
    % background regions untouched by deconvolution recover their real,
    % fluctuating raw value -- a flat +offsetUsed alone would pin every
    % clamped pixel to the same number, i.e. zero background variance.
    sumScaled = sum(decon + negResidual, 3) * opts.scaleFactor + offsetUsed;
else
    sumScaled = sum(decon, 3) * opts.scaleFactor;
end

%% ---- 6b. NON-LINEAR all-in-focus image for TRACING ONLY --------------
% Per-pixel sharpness-weighted mean across z. Weights sum to 1 so values sit in
% single-slice range, but they vary spatially: NOT proportional to signal.
lap = [0 -1 0; -1 4 -1; 0 -1 0];
fm = zeros(size(decon));
for k = 1:nZ
    fm(:,:,k) = imgaussfilt(abs(imfilter(decon(:,:,k), lap, 'replicate')), opts.edfSmooth);
end
w = fm + eps;  w = w ./ sum(w, 3);
if opts.restoreOffset
    segImg = sum((decon + negResidual) .* w, 3) + offsetUsed;
else
    segImg = sum(decon .* w, 3);
end
clear fm w

%% ---- 7. Valid interior (deconvolution ringing margin) ----------------
margin = max(8, 2*rXY);
mask = false(size(sumScaled));
mask(margin+1:end-margin, margin+1:end-margin) = true;

%% ---- 8. uint16 with explicit clipping report -------------------------
% uint16 (not 16-bit float): half-precision steps grow with magnitude (~8 ADU
% near 12,000), i.e. coarsest exactly at filament cores. uint16 gives a uniform
% 1 ADU step and matches the camera's own format.
[R.sum, nClip] = toU16(sumScaled, 'SUM', opts.scaleFactor);
R.seg          = toU16(segImg,    'SEG', 1);
R.mask         = mask;

%% ---- 9. Diagnostics --------------------------------------------------
focus = squeeze(sum(sum(abs(imfilter(img, lap, 'replicate')),1),2));
[~, focusPeakZ] = max(focus);

R.info = struct('channel', string(CHANNEL), 'nZ', nZ, ...
    'offsetUsed', offsetUsed, 'numIter', opts.numIter, ...
    'scaleFactor', opts.scaleFactor, 'doDeconvolution', opts.doDeconvolution, ...
    'fluxRatio', fluxRatio, 'focusPeakZ', focusPeakZ, ...
    'rawMin', rawMin, 'rawMax', rawMax, 'sumMax', double(max(R.sum(:))), ...
    'nSaturatedRaw', nSat, 'nClampedNeg', nClamped, 'nClippedOut', nClip, ...
    'margin', margin, 'segMethod', string(SEG_LABEL), ...
    'pxPerFWHM', pxPerFWHM, 'slicesPerFWHM', szPerFWHM, ...
    'psfSize', string(mat2str(size(psf))), 'firstFile', string(fileList{1}));

%% ---- 10. Optional save ----------------------------------------------
if ~isempty(opts.outBase)
    imwrite(R.sum, [opts.outBase '_SUM_linear.tif']);   % quantify from this
    imwrite(R.seg, [opts.outBase '_SEG.tif']);          % tracing only
    imwrite(uint8(mask)*255, [opts.outBase '_validmask.tif']);
    lo = prctile(double(R.sum(mask)), 1);
    hi = prctile(double(R.sum(mask)), 99.8);
    imwrite(im2uint8(mat2gray(double(R.sum), [lo hi])), [opts.outBase '_display.png']);
    info = R.info; psfUsed = psf; %#ok<NASGU>
    save([opts.outBase '_params.mat'], 'info', 'psfUsed', 'opts', '-v7.3');
end
end

%% ======================= local functions ==============================
function s = setdef(s, f, v)
    if ~isfield(s, f) || isempty(s.(f)), s.(f) = v; end
end

function psf = widefieldPSF(lam_nm, NA, n, dxy, dz, nzPlanes, os)
% Scalar angular-spectrum widefield PSF:
%   PSF(z) = |IFFT{ P(f) .* exp(2i*pi*z*kz(f)) }|^2
% P = circular pupil, kz = axial spatial frequency. Built on an os-times
% oversampled grid then binned, so pixel integration is included - which
% matters because the optical PSF is narrower than one pixel here. Unlike a
% Gaussian this reproduces the defocus cone that actually causes the haze.
%
% NOTE: this is duplicated verbatim in deconvolve_dapi_stack.m. If you ever
% change one, change both, or your two channels stop being comparable.
    lam   = lam_nm/1000;                                            % micron
    coneR = dz*(nzPlanes-1)/2 * tan(asin(min(NA/n,0.999))) / dxy;   % px
    rXY   = ceil(coneR) + 3;                    % PSF must contain the cone
    nb    = 2*rXY + 16;  nb = nb + mod(nb,2);   % even -> centre = nb/2+1
    N     = nb*os;
    d     = dxy/os;

    fx = [0:floor(N/2)-1, -ceil(N/2):-1] / (N*d);
    [FX, FY] = ndgrid(fx, fx);
    FR2   = FX.^2 + FY.^2;
    pupil = double(FR2 <= (NA/lam)^2);
    kz    = sqrt(max((n/lam)^2 - FR2, 0));

    zs    = ((0:nzPlanes-1) - floor(nzPlanes/2)) * dz;
    shift = floor(os/2);                        % os odd -> peak centred in bin
    P     = zeros(nb, nb, nzPlanes);
    for i = 1:nzPlanes
        amp = ifft2(pupil .* exp(2i*pi*zs(i)*kz));
        I   = abs(fftshift(amp)).^2;
        I   = circshift(I, [shift shift]);
        P(:,:,i) = squeeze(sum(sum(reshape(I, os, nb, os, nb), 1), 3));
    end
    c   = nb/2 + 1;
    rXY = min(rXY, c-1);
    psf = P(c-rXY:c+rXY, c-rXY:c+rXY, :);
    psf = psf / sum(psf(:));                    % unit sum -> conserves flux
end

function [out, nOver] = toU16(A, label, sf)
% Round to uint16, reporting clipping rather than hiding it. Silent saturation
% is the easiest way to corrupt a photometric result.
    nOver = sum(A(:) > 65535);
    if nOver > 0
        warning('deconvolve_txred_stack:clip', ...
            '%s: %d px (%.4f%%) exceed 65535 and WILL be clipped. Use scaleFactor <= %.6g (and the SAME value everywhere).', ...
            label, nOver, 100*nOver/numel(A), sf*65535/max(A(:)));
    end
    out = uint16(round(min(max(A,0), 65535)));
end
