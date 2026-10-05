% BATCH_ASMAFILAMENTMASK
% alpha-SMA version of batch_actinFilamentMask (DPP4_20260917 dataset).
% Scans a folder for TxRed (alpha-SMA) SUM projections named
% <Well>_S<##>_TxRed_SUM_linear.tif, segments alpha-SMA stress fibers in each
% with actinFilamentMask, and writes:
%   aSMA_mask_tif\ <Well>_S<##>_TxRed_aSMAmask.tif   uncompressed binary mask
%   aSMA_mask_mat\ <Well>_S<##>_TxRed_aSMAmask.mat   logical mask + metadata
%   aSMA_overlays\ <Well>_S<##>_TxRed_aSMAoverlay.png mask over the alpha-SMA image
%   aSMA_mask_summary.csv                             one row per image
%
% Files are discovered, not assumed — any number of wells/sites is fine, and
% missing sites (e.g. out-of-focus images that were deleted) are simply absent
% from the results rather than causing errors.
%
% Valid-region masks (the deconvolution border exclusion) are picked up
% automatically if present, per site or as one shared file — see below.
%
% actinFilamentMask.m must be on the MATLAB path (or in this folder).

clear; clc;
gpnames={'DDx','DDT','DAT','SDx','SDT','SAT'};    %set

%% ---- settings ----------------------------------------------------------
for ggg=1:size(gpnames,2)
    gpname=gpnames{ggg};
    %dataDir=['D:\Documents\NTU_dropbox\Dropbox\NTU_Pharm\Projects\2026\LabData\DPP4_20260917\',gpname,'/TIFF_files/processed/'];
    dataDir=['E:\浩熏\20260922 HS68(3000) (in SF DMEM + 0.1% BSA) a-SMA Sita(40) TGFb(10) AZD4547 4 days\nd\',gpname,'/TIFF_files/processed/'];
    tifDir  = fullfile(dataDir, 'aSMA_mask_tif');
    matDir  = fullfile(dataDir, 'aSMA_mask_mat');
    pngDir  = fullfile(dataDir, 'aSMA_overlays');
    csvPath = fullfile(dataDir, 'aSMA_mask_summary.csv');

    pixelSize = 0.3642;        % µm/px — from .nd2 metadata (same value used for deconvolution)

    % Valid-region mask. Per-site file <Well>_S<##>_TxRed_validmask.tif is used if
    % it exists; otherwise sharedValidMask is used for every image; otherwise the
    % whole field is analysed. Set sharedValidMask = '' to disable the fallback.
    % The deconvolution script writes a per-site validmask for every stack, so no
    % shared fallback is needed here.
    sharedValidMask = '';

    % What to write
    saveTIF = true;            % uncompressed 1-bit binary TIFF
    saveMAT = true;            % logical mask in a .mat
    savePNG = true;            % RGB overlay for visual QC

    % Overlay appearance
    overlayColor   = [1 0 0];     % red
    overlayDisplay = [1 99.5];    % percentile stretch for the grayscale backdrop

    % Filename contract: <Well>_S<##>_TxRed_SUM_linear.tif
    namePattern = '^(?<well>[A-P]\d{2})_S(?<site>\d+)_TxRed_SUM_linear\.tif$';

    % Options forwarded to actinFilamentMask. Keep them here so every image in the
    % batch is segmented identically — per-field tuning would make the coverage
    % numbers incomparable across conditions. These are the actinFilamentMask
    % defaults, spelled out explicitly so a future default change in the function
    % can't silently change what this batch run used.
    opts = {'PixelSize',      pixelSize, ...
            'CellThreshold',  2, ...          % was 10
            'BackgroundPercentile', 2, ...
            'FilamentWidth',  [0.3 2.5], ...
            'NormWindow',     10, ...
            'StructSens',     2, ...
            'LocalSdCap',     2, ...
            'RidgeThreshold', 0.06, ...       % was 0.04
            'ZFloor',         2.5, ...        % was 1.5
            'MinLength',      8/pixelSize};   % was 3 (µm) — now 8 µm, in px as the fn expects...

    %% ---- discover files ----------------------------------------------------
    if ~isfolder(dataDir)
        warning('Data folder not found - group skipped:\n%s', dataDir);
        continue
    end

    listing = dir(fullfile(dataDir, '*.tif'));
    listing = listing(~[listing.isdir]);

    tok  = regexp({listing.name}, namePattern, 'names', 'once');
    keep = ~cellfun(@isempty, tok);

    nIgnored = nnz(~keep);
    listing  = listing(keep);
    tok      = [tok{keep}];

    if isempty(listing)
        warning(['No files matching the expected name format were found in:\n%s\n' ...
            'Expected e.g. H10_S01_TxRed_SUM_linear.tif - group skipped.'], dataDir);
        continue
    end

    % Sort by well, then site, so the CSV comes out in a sensible order
    wells = string({tok.well})';
    sites = str2double({tok.site})';
    [~, order] = sortrows(table(wells, sites), {'wells', 'sites'});
    listing = listing(order);
    wells   = wells(order);
    sites   = sites(order);

    fprintf('=== %s ===\n', gpname);
    fprintf('Found %d image(s) matching the name format', numel(listing));
    if nIgnored > 0
        fprintf(' (%d other .tif file(s) ignored)', nIgnored);
    end
    fprintf('.\n\n');

    if saveTIF && ~isfolder(tifDir), mkdir(tifDir); end
    if saveMAT && ~isfolder(matDir), mkdir(matDir); end
    if savePNG && ~isfolder(pngDir), mkdir(pngDir); end

    rows = [];   % struct array, one entry per successfully processed image

    %% ---- loop --------------------------------------------------------------
    for i = 1:numel(listing)

        name    = listing(i).name;
        imgPath = fullfile(dataDir, name);
        stem    = sprintf('%s_S%02d_TxRed', wells(i), sites(i));

        fprintf('[%d/%d] %s ... ', i, numel(listing), name);

        % --- valid-region mask: per-site file, else shared, else whole field
        vmPath = fullfile(dataDir, sprintf('%s_validmask.tif', stem));
        if ~isfile(vmPath)
            vmPath = sharedValidMask;
        end
        thisOpts = opts;
        vmSource = 'none';
        if ~isempty(vmPath) && isfile(vmPath)
            try
                vm = imread(vmPath) > 0;
                thisOpts = [opts, {'ValidMask', vm}];      %#ok<AGROW>
                [~, vmName, vmExt] = fileparts(vmPath);
                vmSource = [vmName vmExt];
            catch ME
                warning('Could not read valid mask %s: %s', vmPath, ME.message);
            end
        end

        try
            I = imread(imgPath);
            [mask, info] = actinFilamentMask(I, thisOpts{:});
        catch ME
            fprintf('FAILED\n');
            warning('Failed on %s: %s', name, ME.message);
            continue
        end

        % --- 1-bit binary TIFF, no compression
        maskName = sprintf('%s_aSMAmask.tif', stem);
        if saveTIF
            imwrite(mask, fullfile(tifDir, maskName), 'Compression', 'none');
        end

        % --- .mat: logical mask plus enough metadata to know what it came from
        matName = sprintf('%s_aSMAmask.mat', stem);
        if saveMAT
            maskInfo = struct('SourceFile',   name, ...
                'Group',        gpname, ...
                'Stain',        'alpha-SMA', ...
                'Well',         char(wells(i)), ...
                'Site',         sites(i), ...
                'PixelSize_um', pixelSize, ...
                'ValidMask',    vmSource, ...
                'Options',      {opts}, ...
                'Created',      string(datetime('now')));
            save(fullfile(matDir, matName), 'mask', 'maskInfo', '-v7.3');
        end

        % --- overlay PNG for visual QC
        pngName = sprintf('%s_aSMAoverlay.png', stem);
        if savePNG
            Id   = double(I);
            lims = prctile(Id(info.ValidMask), overlayDisplay);
            if lims(2) <= lims(1)
                lims = [min(Id(:)) max(Id(:))];
            end
            base = im2uint8(mat2gray(Id, lims));
            rgb  = imoverlay(base, mask, overlayColor);
            imwrite(rgb, fullfile(pngDir, pngName));
        end

        % --- collect scalar metrics for the CSV
        cc = bwconncomp(mask);

        r.Group             = string(gpname);
        r.File              = string(name);
        r.Well              = wells(i);
        r.Site              = sites(i);
        r.MaskFile          = string(maskName);
        r.ValidMaskFile     = string(vmSource);
        r.PixelSize_um      = pixelSize;
        r.ValidArea_um2     = nnz(info.ValidMask) * pixelSize^2;
        r.CellArea_um2      = nnz(info.CellMask)  * pixelSize^2;
        r.FilamentArea_um2  = nnz(mask)           * pixelSize^2;
        r.Coverage          = info.Coverage;          % filament / valid area
        r.CoverageInCells   = info.CoverageInCells;   % filament / cell area
        r.CellFraction      = info.CellFraction;      % cell / valid area
        r.NumObjects        = cc.NumObjects;
        r.NoiseSigma        = info.NoiseSigma;
        r.BackgroundLevel   = info.BackgroundLevel;   % watch across conditions --
                                                       % should track each site's
                                                       % real background, not be
                                                       % suspiciously uniform
        r.MedianLocalSd     = info.MedianLocalSd;     % watch across conditions —
        r.SdCap             = info.SdCap;             % see note below
        r.SkeletonLength_um = info.SkeletonLength_um;
        r.ThicknessPx       = string(mat2str(info.ThicknessPx));

        rows = [rows; r];  %#ok<AGROW>

        fprintf('%.3f in-cell coverage, %d objects\n', ...
            info.CoverageInCells, cc.NumObjects);
    end

    %% ---- write CSV ---------------------------------------------------------
    if isempty(rows)
        warning('No images were processed successfully in %s - nothing written.', gpname);
        continue
    end

    T = struct2table(rows);
    writetable(T, csvPath);

    nFailed = numel(listing) - height(T);
    fprintf('\nDone. %d of %d image(s) processed', height(T), numel(listing));
    if nFailed > 0
        fprintf(' (%d failed — see warnings above)', nFailed);
    end
    fprintf('.\n');
    if saveTIF, fprintf('Masks (tif): %s\n', tifDir); end
    if saveMAT, fprintf('Masks (mat): %s\n', matDir); end
    if savePNG, fprintf('Overlays:    %s\n', pngDir); end
    fprintf('CSV:         %s\n\n', csvPath);
end

% NOTE on MedianLocalSd / SdCap: the local-SD cap that rescues bright bundles
% (see actinFilamentMask.m) derives its cap from each field's own median local
% SD. If fiber density differs systematically between conditions (e.g. +/-TGFb), the
% cap — and therefore CoverageInCells — could partly reflect that rather than
% pure biology. Compare MedianLocalSd across groups/conditions before trusting
% CoverageInCells as the headline readout. If it drifts between groups, rerun
% with a fixed 'LocalSdCap' value (or a fixed absolute 'SdCap' override) shared
% across all groups instead of the per-field default.
