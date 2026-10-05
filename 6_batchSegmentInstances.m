%% batchSegmentInstances.m
% Loop over all sites present in a folder and split touching nuclei /
% adjacent cells into labeled instances using segmentInstances.m.
%
% Reads:
%   <inputFolder>/<Well>_S<site>_DAPI_SUM_linear_BgSubtracted.tif
%   <inputFolder>/<Well>_S<site>_TxRed_SUM_linear_BgSubtracted.tif
%   <inputFolder>/binary_masks/<Well>_S<site>_nucleusBinary.tif
%   <inputFolder>/binary_masks/<Well>_S<site>_cellBinary.tif
% (the same folder layout generateBinaryMasks / batchGenerateBinaryMasks
% already produced). Both DAPI and TxRed must be the background-flattened
% '..._BgSubtracted.tif' variant -- same requirement as generateBinaryMasks.m.
%
% Does NOT assume a fixed number of sites -- discovers whatever sites have
% a TxRed image, a matching DAPI image, AND both binary masks present;
% anything missing is skipped with a warning, not an error.
%
% Edit gpnames / groupRange / baseFolder below, then run.

clear; clc;

% ---- USER SETTINGS ----------------------------------------------------
gpnames={'DDx','DDT','DAT','SDx','SDT','SAT'};    %set
groupRange = 1:numel(gpnames); % which groups to process this run, e.g. 1:2 for the a-SMA groups only

% Base project folder, per platform (edit the path for your machine(s)
% below -- the right one is picked automatically via ispc/ismac so you
% don't have to comment/uncomment lines when switching computers).
if ispc
    baseFolder = 'E:\浩熏\20260922 HS68(3000) (in SF DMEM + 0.1% BSA) a-SMA Sita(40) TGFb(10) AZD4547 4 days\nd\';
elseif ismac
    baseFolder = '/Users/FC/Library/CloudStorage/Dropbox/NTU_Pharm/Projects/2026/LabData/DPP4_20260917';
else
    error('batchSegmentInstances:unknownPlatform', ...
        'Add a baseFolder for this platform before running.');
end

%%
for ggg = groupRange
    gpname = gpnames{ggg};
    fprintf('\n=== Group %s (%d/%d) ===\n', gpname, ggg, numel(gpnames));

    % ---- USER SETTINGS ----------------------------------------------------
    inputFolder  = fullfile(baseFolder, gpname, 'TIFF_files', 'processed');
    maskFolder   = fullfile(inputFolder, 'binary_masks');     % contains *_nucleusBinary.tif / *_cellBinary.tif
    outputFolder = fullfile(inputFolder, 'instance_masks');   % labeled masks + QC go here

    txredPattern = '^(?<well>.+)_S(?<site>\d+)_TxRed_SUM_linear_BgSubtracted\.tif$';
    % -------------------------------------------------------------------

    txredFiles = dir(fullfile(inputFolder, '*_TxRed_SUM_linear_BgSubtracted.tif'));
    if isempty(txredFiles)
        warning('No *_TxRed_SUM_linear_BgSubtracted.tif files found in %s -- skipping group.', inputFolder);
        continue;
    end

    siteKeys = {};
    siteInfo = containers.Map('KeyType', 'char', 'ValueType', 'any');
    for i = 1:numel(txredFiles)
        tok = regexp(txredFiles(i).name, txredPattern, 'names');
        if isempty(tok)
            continue;
        end
        key = sprintf('%s_S%s', tok.well, tok.site);
        siteInfo(key) = struct('txred', fullfile(inputFolder, txredFiles(i).name), ...
            'dapi',          fullfile(inputFolder, [key '_DAPI_SUM_linear_BgSubtracted.tif']), ...
            'nucleusBinary', fullfile(maskFolder, [key '_nucleusBinary.tif']), ...
            'cellBinary',    fullfile(maskFolder, [key '_cellBinary.tif']));
        siteKeys{end+1} = key; %#ok<AGROW>
    end
    siteKeys = sort(siteKeys);

    fprintf('Found %d candidate site(s) in %s\n', numel(siteKeys), inputFolder);

    summary = table('Size', [0 4], 'VariableTypes', {'string', 'double', 'double', 'string'}, ...
        'VariableNames', {'Site', 'NumCells', 'NumCellsNoNucleus', 'Status'});

    for i = 1:numel(siteKeys)
        key = siteKeys{i};
        entry = siteInfo(key);

        if ~isfile(entry.nucleusBinary) || ~isfile(entry.cellBinary) || ~isfile(entry.dapi)
            warning('Skipping %s: missing DAPI image or binary mask(s) in %s', key, maskFolder);
            summary = [summary; {key, NaN, NaN, "skipped (missing DAPI or binary mask)"}]; %#ok<AGROW>
            continue;
        end

        fprintf('Processing %s ...\n', key);
        try
            result = segmentInstances(entry.nucleusBinary, entry.cellBinary, entry.dapi, entry.txred, ...
                'OutputDir', outputFolder, ...
                'BaseName', key);
            fprintf('  -> %d cells segmented (%d without a nucleus).\n', ...
                result.numCells, result.numCellsWithoutNucleus);
            summary = [summary; {key, result.numCells, result.numCellsWithoutNucleus, "ok"}]; %#ok<AGROW>
        catch ME
            warning('Failed on %s: %s', key, ME.message);
            fprintf(2, '--- Full error report for %s ---\n%s\n---\n', ...
                key, getReport(ME, 'extended', 'hyperlinks', 'off'));
            summary = [summary; {key, NaN, NaN, "error: " + string(ME.message)}]; %#ok<AGROW>
        end
    end

    disp(summary);

    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end
    writetable(summary, fullfile(outputFolder, 'instance_segmentation_summary.csv'));
    fprintf('Done. Labeled instance masks and QC overlays saved to %s\n', outputFolder);
end
