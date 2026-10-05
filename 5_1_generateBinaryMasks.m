function result = generateBinaryMasks(dapiFile, txredFile, varargin)
%GENERATEBINARYMASKS  Binary (0/1) nucleus and cell foreground masks.
%
%   result = GENERATEBINARYMASKS(dapiFile, txredFile) produces two plain
%   binary masks -- no instance labels, no watershed splitting, no
%   nucleus/cell pairing:
%     nucleusMask : 1 where DAPI signal is real nuclear chromatin
%     cellMask    : 1 where TxRed signal is real cell body (cytoplasm +
%                   nucleus footprint), 0 = true background
%   This intentionally stops BEFORE splitting touching nuclei or
%   separating touching cells, so foreground detection can be checked and
%   tuned on its own, independent of instance-segmentation issues.
%
%   IMPORTANT: BOTH dapiFile and txredFile are now expected to already be
%   regionally background-subtracted (flattened) upstream of this
%   function -- e.g. the '..._BgSubtracted.tif' files from your own
%   regional correction step. Earlier versions of this function tried to
%   handle uncorrected, regionally-varying background internally (first a
%   single global threshold, then a block-wise local threshold map, then
%   a gradient-anchored watershed for TxRed; a rolling-ball imopen
%   subtraction for DAPI). All of that complexity existed only to work
%   around background/glow level genuinely varying across the field --
%   once that's corrected upstream, a single, simple global threshold
%   works reliably again for both channels, and re-running a rolling-ball
%   subtraction on an already-flattened DAPI image would just re-subtract
%   background that's already gone (risking clipping real dim signal).
%   Verified directly for TxRed: on 9 real sites, coverage was stable at
%   56-92% with no degenerate cases, versus wild swings and site-specific
%   artifacts (spurious holes, false positives) with the internal
%   per-region approaches. Feeding this function non-flattened DAPI or
%   TxRed images will likely reproduce the old problems, since the simple
%   thresholds below have no per-region adjustment of their own.
%
%   dapiFile/txredFile may be file paths (char/string) or already-loaded
%   image matrices.
%
%   *** 2026-09-27 UPDATE: simplified nucleus detector ***
%   For clean, high-contrast DAPI images (near-flat background, no real
%   halo/glow to trim -- e.g. dataset a-SMA_20260821), the old adaptive-
%   threshold + per-blob dynamic core-threshold pipeline below was found to
%   FRAGMENT each real nucleus into many small irregular pieces instead of
%   clean ovals: verified directly on site A05_S01, it produced 145 jagged
%   fragments where only ~32 real, clearly separated oval nuclei exist by
%   eye. Root cause: the adaptive-threshold step's local background
%   comparison latches onto faint fibrous background haze across the whole
%   field, merging it with real nuclei into one huge "candidate" blob, and
%   the per-blob core-threshold then carves that blob into many disjoint
%   bright fragments rather than whole nuclei.
%   segmentNucleiMask below now uses a single global Otsu threshold
%   (graythresh) on the smoothed image instead, which is enough when nuclei
%   are much brighter than background and there is no meaningful halo to
%   trim. Verified on A05_S01: imgaussfilt(sigma=2) -> graythresh ->
%   imfill('holes') -> imopen(disk(1)) -> bwareaopen(MinNucleusArea)
%   reproduces exactly 32 nuclei matching the visual count, with clean oval
%   shapes (eccentricity 0.44-0.91, solidity 0.75-0.99) and a clear area gap
%   between debris (<=196 px^2) and real nuclei (>=318 px^2) in that test.
%   For that dataset, also pass 'MinNucleusArea', 300 (up from the default
%   150) -- see the area-gap note above.
%   The OLD adaptive + per-blob core-threshold code is kept below, commented
%   out (prefixed with %), in case a future, noisier/haloed dataset needs it
%   again. 'NucleusSensitivity', 'NucleusCoreFraction',
%   'NucleusCoreSmoothSigma', 'NucleusPeakPercentile', and
%   'NucleusErosionRadius' below are now UNUSED by the active code path --
%   they only apply if you restore the commented-out block.
%
%   Name-Value options:
%     'SmoothSigma'             gaussian smoothing sigma (px). Default 2.
%     'NucleusSensitivity'      [LEGACY -- unused by the active simplified
%                               detector, see update note above] adaptive
%                               threshold sensitivity for DAPI, 0-1, higher
%                               = more permissive. Default 0.45.
%     'MinNucleusArea'          minimum connected-component area (px^2)
%                               to keep as nucleus signal. Default 150.
%                               For clean/high-contrast DAPI (see update
%                               note above), pass 300 instead.
%     'NucleusCoreFraction'     [LEGACY -- unused by the active simplified
%                               detector, see update note above]
%                               the adaptive threshold above is only used
%                               to LOCATE candidate nucleus blobs -- it is
%                               a local-contrast criterion, so it also
%                               picks up a dim halo of real-but-faint
%                               signal well beyond each nucleus's bright
%                               core, and that halo's width varies a lot
%                               from nucleus to nucleus (checked directly:
%                               one real nucleus's signal didn't fall to
%                               background until >30px out from a core
%                               that was only ~15px). A single fixed-size
%                               erosion can't track that variable halo
%                               width -- correctly trimming a wide-halo
%                               nucleus under-trims a narrow-halo one, or
%                               vice versa. Instead, EACH candidate blob is
%                               re-thresholded using its own dynamic range:
%                               keep only pixels at or above
%                               dapiBackgroundLevel + NucleusCoreFraction x
%                               (blobPeakIntensity - dapiBackgroundLevel),
%                               where blobPeakIntensity is that blob's own
%                               99th-percentile raw intensity. 0-1, higher
%                               = tighter/smaller mask. Default 0.5. Tune
%                               this first if masks still look too
%                               big/small -- it has far more effect than
%                               'NucleusErosionRadius' below.
%     'NucleusCoreSmoothSigma'  [LEGACY -- unused by the active simplified
%                               detector, see update note above]
%                               gaussian sigma (px) applied to the DAPI image
%                               used ONLY for the per-blob core threshold
%                               above (the blob peak and the pixel-wise
%                               comparison). Default 0 (raw pixels, the
%                               original behaviour). Deconvolved DAPI has
%                               sharp chromatin texture: comparing raw
%                               pixels against a peak set by the brightest
%                               speckles keeps only the speckles, so each
%                               nucleus breaks into small jagged fragments
%                               (several seeds per nucleus, or none once
%                               fragments fall under MinNucleusArea). ~2
%                               restores whole, smooth nuclei on
%                               deconvolved data.
%     'NucleusPeakPercentile'   [LEGACY -- unused by the active simplified
%                               detector, see update note above]
%                               percentile of the blob's (optionally
%                               smoothed) intensities used as its peak for
%                               the core threshold. Default 99.
%     'NucleusErosionRadius'    [LEGACY -- unused by the active simplified
%                               detector, which uses a small imopen(disk(1))
%                               instead; see update note above]
%                               disk SE radius (px) for a small final
%                               erosion, purely to smooth the jagged pixel
%                               boundary the per-blob core threshold above
%                               leaves behind -- not meant to do the main
%                               shrinking (that's 'NucleusCoreFraction').
%                               Default 1.
%     'MinNucleusMeanIntensity' minimum mean raw DAPI intensity, ABOVE the
%                               estimated DAPI background level (see
%                               'BackgroundPercentile'), a candidate
%                               nucleus blob must have to be kept (rejects
%                               faint debris/dust that passes the adaptive
%                               threshold on local contrast alone). Default
%                               50. This MUST be measured relative to
%                               background, not as an absolute intensity:
%                               once upstream processing restores a real
%                               camera baseline (e.g. ~400-600 ADU rather
%                               than ~0), an absolute 50 is below the
%                               baseline itself and rejects nothing.
%     'BackgroundPercentile'    percentile (0-100) of pixel intensities
%                               assumed to be true background -- used for
%                               both the DAPI background level (debris
%                               rejection above) and the (now-flattened)
%                               TxRed background level/noise for the cell
%                               threshold. Default 2 (tightened from an
%                               initial 10: the "background" pool in a
%                               dense field gets contaminated by real dim
%                               cytoplasm at looser percentiles).
%     'CellThresholdNumSigma'   cell foreground = pixels more than this
%                               many background-noise standard deviations
%                               above the TxRed background level. Default
%                               8. Verified against 9 real (background-
%                               flattened) sites: gives coverage in a sane
%                               56-92% range everywhere, and correctly
%                               excludes every false-positive point
%                               reported against earlier versions, with
%                               comfortable margin (not borderline).
%     'MinCellMarginAbsolute'   floor on the margin above background,
%                               regardless of how small CellThresholdNumSigma
%                               x bgNoise comes out to. Default 150.
%                               Needed because bgNoise itself can vary a
%                               lot between sites even after background
%                               flattening -- one real site measured
%                               bgNoise=4.8 (vs 24-41 on 8 others), so its
%                               sigma-based margin was only ~39 while
%                               others sat at 190-326. That tiny margin let
%                               a dim, roughly-constant-width out-of-focus
%                               halo around each cell get counted as cell,
%                               pushing the mask boundary outside the real
%                               edge on that one unusually clean site. The
%                               floor keeps outlier-low-noise sites from
%                               getting an outlier-low threshold.
%     'MinCellArea'             minimum connected-component area (px^2)
%                               to keep as cell signal. Default 400.
%     'MaxHoleFillArea'         largest enclosed background hole (px^2)
%                               that will be filled in as cell signal.
%                               Larger enclosed gaps (fully surrounded by
%                               cells but not touching the image border)
%                               are left as background, since they can be
%                               genuine extracellular gaps rather than a
%                               dim patch inside one cell. Default 300
%                               (tightened from an initial 500 after a
%                               ~660 px^2 hole was still getting filled).
%     'OutputDir'               if provided, writes
%                               '<BaseName>_nucleusBinary.tif' and
%                               '<BaseName>_cellBinary.tif' (uint8, values
%                               0/1) plus a QC overlay PNG into this
%                               folder. Default '' (nothing written).
%     'BaseName'                base filename for saved outputs. Default
%                               derived from dapiFile if it is a path,
%                               otherwise 'site'.
%     'SaveQCOverlay'           true/false, save a mask-overlay PNG when
%                               OutputDir is set. Default true.
%
%   Output struct fields:
%     .nucleusMask   logical, 1608x1608 (or input size)
%     .cellMask      logical
%     .dapiFile / .txredFile   echoed inputs
%     .dapiBackgroundLevel               estimated DAPI background level
%     .txredBackgroundLevel, .txredBackgroundNoise
%                    estimated TxRed background level/noise -- watch
%                    these across a batch; they should now be far more
%                    consistent site-to-site than before your regional
%                    background correction was added upstream

p = inputParser;
addParameter(p, 'SmoothSigma', 2);
addParameter(p, 'NucleusSensitivity', 0.45);
addParameter(p, 'MinNucleusArea', 150);
addParameter(p, 'NucleusCoreFraction', 0.5);
addParameter(p, 'NucleusCoreSmoothSigma', 0);
addParameter(p, 'NucleusPeakPercentile', 99);
addParameter(p, 'NucleusErosionRadius', 1);
addParameter(p, 'MinNucleusMeanIntensity', 50);
addParameter(p, 'BackgroundPercentile', 2);
addParameter(p, 'CellThresholdNumSigma', 8);
addParameter(p, 'MinCellMarginAbsolute', 150);
addParameter(p, 'MinCellArea', 400);
addParameter(p, 'MaxHoleFillArea', 300);
addParameter(p, 'OutputDir', '');
addParameter(p, 'BaseName', '');
addParameter(p, 'SaveQCOverlay', true);
parse(p, varargin{:});
opt = p.Results;

[dapiImg, dapiPathStr] = loadImage(dapiFile);
[txredImg, ~]          = loadImage(txredFile);

if ~isequal(size(dapiImg), size(txredImg))
    error('generateBinaryMasks:sizeMismatch', ...
        'DAPI and TxRed images must be the same size (got %s vs %s).', ...
        mat2str(size(dapiImg)), mat2str(size(txredImg)));
end

if isempty(opt.BaseName)
    if ~isempty(dapiPathStr)
        [~, baseName] = fileparts(dapiPathStr);
        baseName = regexprep(baseName, '_DAPI.*$', '');
    else
        baseName = 'site';
    end
else
    baseName = opt.BaseName;
end

% ---- nucleus binary mask from DAPI ------------------------------------
dapiBackgroundLevel = estimateBackgroundLevel(dapiImg, opt);
nucleusBW = segmentNucleiMask(dapiImg, dapiBackgroundLevel, opt);
nucleusBW = rejectDebrisByIntensity(nucleusBW, dapiImg, dapiBackgroundLevel, opt);

% ---- cell binary mask from (already-flattened) TxRed ------------------
[cellBW, txredBackgroundLevel, txredBackgroundNoise] = segmentCellForeground(txredImg, nucleusBW, opt);

result.nucleusMask = nucleusBW;
result.cellMask    = cellBW;
result.dapiFile     = dapiFile;
result.txredFile    = txredFile;
result.dapiBackgroundLevel   = dapiBackgroundLevel;
result.txredBackgroundLevel  = txredBackgroundLevel;
result.txredBackgroundNoise  = txredBackgroundNoise;

% ---- optional outputs ---------------------------------------------
if ~isempty(opt.OutputDir)
    if ~exist(opt.OutputDir, 'dir')
        mkdir(opt.OutputDir);
    end
    imwrite(uint8(nucleusBW), fullfile(opt.OutputDir, [baseName '_nucleusBinary.tif']));
    imwrite(uint8(cellBW),    fullfile(opt.OutputDir, [baseName '_cellBinary.tif']));

    if opt.SaveQCOverlay
        saveBinaryQCOverlay(txredImg, dapiImg, cellBW, nucleusBW, ...
            fullfile(opt.OutputDir, [baseName '_binaryQC_overlay.png']));
    end
end

end % generateBinaryMasks


% ======================================================================
% helper functions
% ======================================================================

function [img, pathStr] = loadImage(input)
pathStr = '';
if ischar(input) || isstring(input)
    pathStr = char(input);
    img = imread(pathStr);
else
    img = input;
end
img = double(img);
end


function nucleusBW = segmentNucleiMask(dapiImg, dapiBackgroundLevel, opt) %#ok<INUSL>
% dapiImg is expected to already be regionally background-flattened
% upstream (same assumption as the TxRed path) -- so no internal
% rolling-ball subtraction is done here; doing one on top of an
% already-flattened image would just re-subtract background that's
% already gone and risk clipping real dim nuclear signal.
% Note: this does NOT split touching nuclei -- a clump of several
% touching nuclei stays as one connected foreground blob here.
%
% *** 2026-09-27: simplified global-threshold detector (see the UPDATE
% note in the function-level help above for why). dapiBackgroundLevel is
% no longer used by this active path -- kept as an input argument only so
% the OLD code below can be restored without changing the call site.
imgNorm = mat2gray(dapiImg);
imgSmooth = imgaussfilt(imgNorm, opt.SmoothSigma);

nucleusBW = imbinarize(imgSmooth, graythresh(imgSmooth));
nucleusBW = imfill(nucleusBW, 'holes');
nucleusBW = imopen(nucleusBW, strel('disk', 1));   % just smooths jagged pixel edges
nucleusBW = bwareaopen(nucleusBW, opt.MinNucleusArea);

% ======================================================================
% OLD (adaptive threshold + per-blob dynamic core threshold) -- kept for
% reference / in case a future noisier or haloed dataset needs it again.
% See the UPDATE note in the function-level help above for why this was
% replaced: on clean, high-contrast DAPI it fragmented each real nucleus
% into many small irregular pieces instead of whole ovals.
% ======================================================================
% % Smooth, adaptively threshold, and clean up to LOCATE candidate nucleus
% % blobs -- this stage is a local-contrast criterion, so it also captures
% % a dim halo of real-but-faint signal well beyond each nucleus's bright
% % core (halo width varies a lot nucleus-to-nucleus, so it can't be fixed
% % with a single flat erosion).
% imgNorm = mat2gray(dapiImg);
% imgSmooth = imgaussfilt(imgNorm, opt.SmoothSigma);
%
% candidateBW = imbinarize(imgSmooth, 'adaptive', ...
%     'ForegroundPolarity', 'bright', 'Sensitivity', opt.NucleusSensitivity);
%
% candidateBW = imfill(candidateBW, 'holes');
% candidateBW = imopen(candidateBW, strel('disk', 2));
% candidateBW = bwareaopen(candidateBW, opt.MinNucleusArea);
%
% % Re-threshold each candidate blob using its OWN dynamic range: keep only
% % pixels at or above dapiBackgroundLevel + NucleusCoreFraction x (that
% % blob's own peak - dapiBackgroundLevel). This shrinks each nucleus by
% % however much its own halo actually extends, rather than by a fixed pixel
% % count that over-trims some nuclei and under-trims others.
% % Optionally do this comparison on a smoothed copy (see
% % 'NucleusCoreSmoothSigma') so chromatin texture in deconvolved data doesn't
% % fragment each nucleus.
% if opt.NucleusCoreSmoothSigma > 0
%     coreImg = imgaussfilt(dapiImg, opt.NucleusCoreSmoothSigma);
% else
%     coreImg = dapiImg;
% end
% cc = bwconncomp(candidateBW);
% nucleusBW = false(size(candidateBW));
% for i = 1:cc.NumObjects
%     idx = cc.PixelIdxList{i};
%     blobPeak = prctile(coreImg(idx), opt.NucleusPeakPercentile);
%     coreThreshold = dapiBackgroundLevel + ...
%         opt.NucleusCoreFraction * (blobPeak - dapiBackgroundLevel);
%     keepIdx = idx(coreImg(idx) >= coreThreshold);
%     nucleusBW(keepIdx) = true;
% end
%
% nucleusBW = imfill(nucleusBW, 'holes');
% nucleusBW = bwareaopen(nucleusBW, max(round(opt.MinNucleusArea / 2), 1));
%
% % Small final erosion purely to smooth the jagged boundary the per-blob
% % intensity threshold above leaves behind -- NucleusCoreFraction does the
% % actual shrinking now, this is just polish.
% if opt.NucleusErosionRadius > 0
%     nucleusBW = imerode(nucleusBW, strel('disk', opt.NucleusErosionRadius));
% end
end


function [bgLevel, bgNoise] = estimateBackgroundLevel(img, opt)
% Robust background level/noise from the darkest fraction of pixels.
lowVal = prctile(img(:), opt.BackgroundPercentile);
bgPixels = img(img <= lowVal);
bgLevel = median(bgPixels);
bgNoise = std(bgPixels);
if bgNoise == 0
    bgNoise = eps; % guard against a perfectly flat background patch
end
end


function bw = rejectDebrisByIntensity(bw, rawImg, backgroundLevel, opt)
% Drop connected components whose mean raw intensity, ABOVE the estimated
% background level, is too low to be real signal (faint debris/dust that
% only passed the adaptive threshold because of relative, not absolute,
% contrast). Comparing above background rather than to an absolute
% constant keeps this meaningful regardless of where the camera's real
% baseline happens to sit.
cc = bwconncomp(bw);
keep = false(size(bw));
for i = 1:cc.NumObjects
    idx = cc.PixelIdxList{i};
    if (mean(rawImg(idx)) - backgroundLevel) >= opt.MinNucleusMeanIntensity
        keep(idx) = true;
    end
end
bw = keep;
end


function [cellFG, bgLevel, bgNoise] = segmentCellForeground(txredImg, nucleusBW, opt)
% Whole-cell foreground mask from an already regionally-flattened TxRed
% image. A single global threshold is enough here BECAUSE the input has
% already had its regional background variation corrected upstream --
% see the function-level help for why this replaced several more
% complicated internal approaches that were working around that same
% variation.
imgSmooth = imgaussfilt(txredImg, opt.SmoothSigma);

[bgLevel, bgNoise] = estimateBackgroundLevel(imgSmooth, opt);
margin = max(opt.CellThresholdNumSigma * bgNoise, opt.MinCellMarginAbsolute);
threshold = bgLevel + margin;

cellFG = imgSmooth > threshold;
cellFG = cellFG | nucleusBW;

cellFG = imclose(cellFG, strel('disk', 5));

% Fill only SMALL enclosed holes (e.g. a dim patch inside one cell's
% footprint) -- do not blindly fill every enclosed background region.
% In a confluent monolayer, a large hole fully surrounded by cells on
% all sides can still be a genuine extracellular gap, not part of any
% cell, even though it never touches the image border. Filling
% everything would wrongly turn such gaps into "cell".
filledFull = imfill(cellFG, 'holes');
holes = filledFull & ~cellFG;
holeCC = bwconncomp(holes);
for i = 1:holeCC.NumObjects
    idx = holeCC.PixelIdxList{i};
    if numel(idx) <= opt.MaxHoleFillArea
        cellFG(idx) = true;
    end
end

cellFG = bwareaopen(cellFG, opt.MinCellArea);
end


function saveBinaryQCOverlay(txredImg, dapiImg, cellBW, nucleusBW, outPath)
% Quick visual check: TxRed shown as grayscale/white, DAPI as magenta,
% with cell foreground boundary (yellow) and nucleus foreground boundary
% (cyan) overlaid on top.
txredNorm = mat2gray(imadjust(mat2gray(txredImg)));
dapiNorm  = mat2gray(imadjust(mat2gray(dapiImg)));

% TxRed contributes equally to R/G/B (grayscale/white). DAPI adds to R
% and B only (magenta), so DAPI-only areas read magenta, TxRed-only
% areas read white/gray, and overlap reads a lighter pink-white.
rC = min(txredNorm + dapiNorm, 1);
gC = txredNorm;
bC = min(txredNorm + dapiNorm, 1);

cellBoundary = bwperim(cellBW);
nucBoundary  = bwperim(nucleusBW);

rC(cellBoundary) = 1; gC(cellBoundary) = 1; bC(cellBoundary) = 0; % yellow
rC(nucBoundary)  = 0; gC(nucBoundary)  = 1; bC(nucBoundary)  = 1; % cyan
rgb = cat(3, rC, gC, bC);

imwrite(rgb, outPath);
end
