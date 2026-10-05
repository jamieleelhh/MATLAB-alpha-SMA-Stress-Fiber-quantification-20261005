%% ND2 -> TIFF batch converter
% Reads every .nd2 file in an input folder and writes out TIFF images.
%
% REQUIREMENT: the Bio-Formats MATLAB toolbox ("bfmatlab") must be on the path.
%   1. Download "bfmatlab.zip" from
%      https://www.openmicroscopy.org/bio-formats/downloads/
%   2. Unzip it somewhere permanent.
%   3. In MATLAB run:  addpath('C:\path\to\bfmatlab')   (adjust to your path)
%
% Bio-Formats reads ND2 pixel data at full bit depth (usually 16-bit), so the
% exported TIFFs preserve the original intensity values.

clear; clc;
gpnames={'DDx','DDT','DAT','SDx','SDT','SAT'};    %set

for ggg=1:size(gpnames,2)
    gpname=gpnames{ggg};
    cwd=['E:\浩熏\20260922 HS68(3000) (in SF DMEM + 0.1% BSA) a-SMA Sita(40) TGFb(10) AZD4547 4 days\nd\',gpname,'\'];    %set
    
    %% ===================== USER SETTINGS =====================
    inputFolder  = cwd;    % folder containing your .nd2 files
    outputFolder = [cwd,'TIFF_files\'];  % folder where TIFFs will be written
    saveMode     = 'planes';   % 'planes' = one TIFF per image plane (safest, general)
                               % 'stack'  = one multi-page TIFF per series
                               %            (best for single-channel Z-stacks/timelapses)
    %% =========================================================
    
    % --- verify Bio-Formats is available ---
    if exist('bfGetReader', 'file') ~= 2
        error('Bio-Formats (bfmatlab) not found on the MATLAB path. See the header for setup.');
    end
    
    if ~exist(outputFolder, 'dir'); mkdir(outputFolder); end
    
    nd2Files = dir(fullfile(inputFolder, '*.nd2'));
    if isempty(nd2Files)
        error('No .nd2 files found in: %s', inputFolder);
    end
    fprintf('Found %d ND2 file(s).\n', numel(nd2Files));
    
    for f = 1:numel(nd2Files)
        nd2Path = fullfile(nd2Files(f).folder, nd2Files(f).name);
        [~, baseName] = fileparts(nd2Files(f).name);
        fprintf('\n[%d/%d] %s\n', f, numel(nd2Files), nd2Files(f).name);
    
        reader  = bfGetReader(nd2Path);
        nSeries = reader.getSeriesCount();
    
        for s = 0:nSeries-1
            reader.setSeries(s);
            sizeZ = reader.getSizeZ();
            sizeC = reader.getSizeC();
            sizeT = reader.getSizeT();
    
            % In 'stack' mode everything for this series goes into one multi-page file
            if strcmpi(saveMode, 'stack')
                stackName  = makeName(baseName, nSeries, s, [], [], [], sizeZ, sizeC, sizeT);
                stackPath  = fullfile(outputFolder, [stackName '.tif']);
                firstWrite = true;
            end
    
            for t = 0:sizeT-1
                for z = 0:sizeZ-1
                    for c = 0:sizeC-1
                        iPlane = reader.getIndex(z, c, t) + 1;   % Bio-Formats is 0-based
                        img = bfGetPlane(reader, iPlane);        % full-bit-depth pixels
    
                        switch lower(saveMode)
                            case 'planes'
                                name = makeName(baseName, nSeries, s, t, z, c, sizeZ, sizeC, sizeT);
                                imwrite(img, fullfile(outputFolder, [name '.tif']), ...
                                        'Compression', 'none');
                            case 'stack'
                                if firstWrite
                                    imwrite(img, stackPath, 'Compression', 'none');
                                    firstWrite = false;
                                else
                                    imwrite(img, stackPath, 'Compression', 'none', ...
                                            'WriteMode', 'append');
                                end
                            otherwise
                                error('saveMode must be ''planes'' or ''stack''.');
                        end
                    end
                end
            end
            fprintf('   series %d/%d: %d plane(s)  [Z=%d C=%d T=%d]\n', ...
                    s+1, nSeries, sizeZ*sizeC*sizeT, sizeZ, sizeC, sizeT);
        end
        reader.close();
    end
    
    fprintf('\nDone. TIFFs are in: %s\n\n', outputFolder);
end

%% ---------- helper: build a tidy file name ----------
% Only appends a dimension tag when that dimension actually has >1 plane,
% so simple single-image ND2 files just produce "<basename>.tif".
function name = makeName(base, nSeries, s, t, z, c, sizeZ, sizeC, sizeT)
    name = base;
    if nSeries > 1;              name = sprintf('%s_S%02d', name, s+1); end
    if ~isempty(t) && sizeT > 1; name = sprintf('%s_T%03d', name, t+1); end
    if ~isempty(z) && sizeZ > 1; name = sprintf('%s_Z%03d', name, z+1); end
    if ~isempty(c) && sizeC > 1; name = sprintf('%s_C%02d', name, c+1); end
end