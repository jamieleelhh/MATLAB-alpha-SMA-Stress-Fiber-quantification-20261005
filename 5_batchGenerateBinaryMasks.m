%% batchGenerateBinaryMasks.m
% Loop over all sites present in a single folder and generate binary
% (0/1) nucleus + cell foreground masks for each one using
% generateBinaryMasks.m.
%
% Does NOT assume a fixed number of sites (e.g. 9) -- some sites may be
% missing because off-focus images were deleted. It discovers whatever
% DAPI/TxRed pairs actually exist by parsing filenames.
%
% Expected filename pattern (per current dataset):
%   <Well>_S<site>_DAPI_SUM_linear_BgSubtracted.tif
%   <Well>_S<site>_TxRed_SUM_linear_BgSubtracted.tif
%   Both channels must now be the regionally background-flattened
%   variant -- plain, unsubtracted DAPI/TxRed files are ignored; see
%   generateBinaryMasks.m for why.
%   e.g. H10_S01_DAPI_SUM_linear_BgSubtracted.tif /
%        H10_S01_TxRed_SUM_linear_BgSubtracted.tif
%
% Edit inputFolder / outputFolder below, then run.

clear; clc;

% ---- USER SETTINGS ----------------------------------------------------
gpnames={'DDx','DDT','DAT','SDx','SDT','SAT'};    %set
groupRange = 1:numel(gpnames); % which groups to process this run, e.g. 1:2 for the a-SMA groups only

% Options forwarded to generateBinaryMasks for every site. Only the
% cell-mask margin floor differs from the function defaults:
%   MinCellMarginAbsolute: default 150 ADU was set for the Huh-7 data. In
%   this dataset the TxRed signal sits only ~20-200 ADU above a ~600 ADU
%   background (bg noise after smoothing ~2 ADU), so a 150 floor kept only
%   the bright perinuclear cores (~17% coverage on H10_S01/S02). 20 ADU
%   follows the real cell edges (~72-80% coverage); 30 already cuts off
%   dim cells. The 8-sigma rule (~13-16 ADU here) is still applied -- the
%   floor only stops unusually clean sites from getting a lower threshold.
%   NucleusCoreSmoothSigma / NucleusPeakPercentile / NucleusCoreFraction:
%   the deconvolved DAPI has sharp chromatin texture, and the default
%   raw-pixel core threshold (50% of each blob's 99th percentile) broke
%   every nucleus into small jagged fragments (cyan squiggles in the QC
%   overlays). Thresholding a 2-px-smoothed copy against its 95th
%   percentile at 40% gives whole nuclei with smooth outlines.
maskOpts = {'MinCellMarginAbsolute', 20, ...
            'CellThresholdNumSigma', 8, ...
            'NucleusCoreSmoothSigma', 2, ...
            'NucleusPeakPercentile', 95, ...
            'NucleusCoreFraction', 0.4};

% Base project folder, per platform (edit the path for your machine(s)
% below -- the right one is picked automatically via ispc/ismac so you
% don't have to comment/uncomment lines when switching computers).
if ispc
    baseFolder = 'E:\浩熏\20260922 HS68(3000) (in SF DMEM + 0.1% BSA) a-SMA Sita(40) TGFb(10) AZD4547 4 days\nd\';
elseif ismac
    baseFolder = '/Users/FC/Library/CloudStorage/Dropbox/NTU_Pharm/Projects/2026/LabData/DPP4_20260917';
else
    error('batchSegmentCells:unknownPlatform', ...
        'Add a baseFolder for this platform before running.');
end

%%
for ggg = groupRange
    gpname = gpnames{ggg};
    fprintf('\n=== Group %s (%d/%d) ===\n', gpname, ggg, numel(gpnames));


    % ---- USER SETTINGS ----------------------------------------------------
    inputFolder  = fullfile(baseFolder, gpname, 'TIFF_files', 'processed');
    outputFolder = fullfile(inputFolder, 'binary_masks'); % masks + QC go here

    filePattern  = '^(?<well>.+)_S(?<site>\d+)_(?<channel>DAPI|TxRed)_SUM_linear(?<bgsub>_BgSubtracted)?\.tif$';
    % Both nucleus- and cell-mask quality now depend on being fed an
    % already regionally background-flattened image (see
    % generateBinaryMasks.m) -- so DAPI and TxRed files are only picked
    % up here if they carry the '_BgSubtracted' suffix; the plain,
    % unflattened files of either channel are ignored.
    % -------------------------------------------------------------------

    allFiles = dir(fullfile(inputFolder, '*_SUM_linear*.tif'));
    if isempty(allFiles)
        warning('No *_SUM_linear*.tif files found in %s - group skipped.', inputFolder);
        continue
    end

    % Group files by well+site, recording which channel(s) were found.
    siteMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
    for i = 1:numel(allFiles)
        tok = regexp(allFiles(i).name, filePattern, 'names');
        if isempty(tok)
            continue; % skip files that don't match the expected naming
        end
        isBgSub = ~isempty(tok.bgsub);
        if ~isBgSub
            continue; % skip unflattened files of either channel -- see note above
        end

        key = sprintf('%s_S%s', tok.well, tok.site);
        if ~isKey(siteMap, key)
            siteMap(key) = struct('well', tok.well, 'site', tok.site, ...
                'dapi', '', 'txred', '');
        end
        entry = siteMap(key);
        if strcmpi(tok.channel, 'DAPI')
            entry.dapi = fullfile(inputFolder, allFiles(i).name);
        else
            entry.txred = fullfile(inputFolder, allFiles(i).name);
        end
        siteMap(key) = entry;
    end

    siteKeys = sort(keys(siteMap));
    fprintf('Found %d candidate site(s) in %s\n', numel(siteKeys), inputFolder);

    summary = table('Size', [0 8], ...
        'VariableTypes', {'string', 'double', 'double', 'double', 'double', 'double', 'double', 'string'}, ...
        'VariableNames', {'Site', 'NucleusFgFraction', 'CellFgFraction', ...
                          'DapiBackgroundLevel', 'TxRedBackgroundLevel', 'TxRedBackgroundNoise', ...
                          'TxRedMarginUsed', 'Status'});
    optStruct = struct(maskOpts{:});

    for i = 1:numel(siteKeys)
        key = siteKeys{i};
        entry = siteMap(key);

        if isempty(entry.dapi) || isempty(entry.txred)
            warning('Skipping %s: missing %s channel.', key, ...
                iif(isempty(entry.dapi), 'DAPI', 'TxRed'));
            summary = [summary; {key, NaN, NaN, NaN, NaN, NaN, NaN, "skipped (missing channel)"}]; %#ok<AGROW>
            continue;
        end

        fprintf('Processing %s ...\n', key);
        try
            result = generateBinaryMasks(entry.dapi, entry.txred, maskOpts{:}, ...
                'OutputDir', outputFolder, ...
                'BaseName', key);
            % margin actually applied (same rule as inside generateBinaryMasks):
            % equals MinCellMarginAbsolute whenever the floor is the binding term
            marginUsed = max(optStruct.CellThresholdNumSigma * result.txredBackgroundNoise, ...
                             optStruct.MinCellMarginAbsolute);
            nucFrac = mean(result.nucleusMask(:));
            cellFrac = mean(result.cellMask(:));
            fprintf('  -> nucleus fg %.1f%%, cell fg %.1f%% (bg: DAPI %.0f, TxRed %.0f +/- %.1f, margin %.1f)\n', ...
                nucFrac*100, cellFrac*100, result.dapiBackgroundLevel, ...
                result.txredBackgroundLevel, result.txredBackgroundNoise, marginUsed);
            summary = [summary; {key, nucFrac, cellFrac, result.dapiBackgroundLevel, ...
                result.txredBackgroundLevel, result.txredBackgroundNoise, marginUsed, "ok"}]; %#ok<AGROW>
        catch ME
            warning('Failed on %s: %s', key, ME.message);
            summary = [summary; {key, NaN, NaN, NaN, NaN, NaN, NaN, "error: " + string(ME.message)}]; %#ok<AGROW>
        end
    end

    disp(summary);

    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end
    writetable(summary, fullfile(outputFolder, 'binary_mask_summary.csv'));
    fprintf('Done. Binary masks and QC overlays saved to %s\n', outputFolder);
end

function out = iif(cond, a, b)
% small inline ternary helper for the warning message above
if cond
    out = a;
else
    out = b;
end
end
