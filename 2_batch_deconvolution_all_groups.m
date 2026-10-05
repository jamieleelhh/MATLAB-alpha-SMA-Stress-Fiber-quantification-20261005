%%
close all;clear;clc
gpnames={'DDx','DDT','DAT','SDx','SDT','SAT'};    %set

%% Acquisition settings for this image series
nZexp = 9;      % slices per z-stack (was 9 in earlier series)
zStep = 2.40;    % micron between slices (confirmed from .nd2 metadata)
pixelSize = 0.3642;  % micron per pixel (confirmed from .nd2 metadata)

%%
for ggg=1:size(gpnames,2)
    gpname=gpnames{ggg};
%    cwd=['D:\Documents\NTU_dropbox\Dropbox\NTU_Pharm\Projects\2026\LabData\DPP4_20260917\',gpname,'/TIFF_files/'];
    cwd=['E:\浩熏\20260922 HS68(3000) (in SF DMEM + 0.1% BSA) a-SMA Sita(40) TGFb(10) AZD4547 4 days\nd\',gpname,'/TIFF_files/'];
    channelnames={'TxRed','DAPI'};offsetval=[601,414];

    %%
    for currentch=1:2

        %%
        opts = struct('NA',0.75,'n_medium',1.00,'pixelSize',pixelSize,'zStep',zStep, ...
            'cameraOffset',offsetval(currentch), 'restoreOffset', true, ...   % <- offset from a DARK FRAME
            'scaleFactor', 1/nZexp);   % z-sum -> z-average; keep constant within this experiment

        inFolder  = cwd;
        outFolder = fullfile(inFolder,'processed');
        if ~exist(outFolder,'dir'), mkdir(outFolder); end

        d   = dir(fullfile(inFolder,['Well*_Channel',channelnames{currentch},'_Seq*_Z*.tif']));
        tok = regexp({d.name}, ...
            '^Well([^_]+)_Channel([^_]+)_Seq(\d+)_S(\d+)_Z(\d+)\.tif$','tokens','once');
        keep = ~cellfun(@isempty, tok);   d = d(keep);   tok = tok(keep);
        key  = cellfun(@(t) sprintf('%s_S%s', t{1}, t{4}), tok, 'UniformOutput', false);
        seq  = cellfun(@(t) t{3}, tok, 'UniformOutput', false);

        uk = unique(key);   rows = [];
        for i = 1:numel(uk)
            disp([uk{i},'_',channelnames{currentch}])
            sel = strcmp(key, uk{i});
            f   = fullfile({d(sel).folder}, {d(sel).name});
            nSeq = numel(unique(seq(sel)));
            if nSeq > 1                                   % field acquired more than once
                warning('%s: acquired in %d Seq - skipped.', uk{i}, nSeq);
                continue
            end
            if numel(f) ~= nZexp
                warning('%s: %d slices, expected %d - skipped.', uk{i}, numel(f), nZexp);
                continue
            end
            opts.outBase = fullfile(outFolder, [uk{i} '_',channelnames{currentch}]);
            if currentch==1
                R = deconvolve_txred_stack(f, opts);
            elseif currentch==2
                R = deconvolve_dapi_stack(f, opts);
            end
            rows = [rows; R.info];
        end
        if isempty(rows)
            warning('No stacks processed for %s - no manifest written.', channelnames{currentch});
        else
            writetable(struct2table(rows), fullfile(outFolder,[channelnames{currentch},'_manifest.csv']));
        end
    end
end
