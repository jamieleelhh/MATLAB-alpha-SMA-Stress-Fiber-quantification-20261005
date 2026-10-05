%%
clear;clc;close all

%% cwd -- pick automatically based on OS so the same script runs on Mac and PC
if ispc
    cwd='E:\浩熏\20260922 HS68(3000) (in SF DMEM + 0.1% BSA) a-SMA Sita(40) TGFb(10) AZD4547 4 days\nd\';
else
    cwd='/Users/FC/Library/CloudStorage/Dropbox/NTU_Pharm/Projects/2026/LabData/DPP4_20260917/';
end
% strip any trailing separator so fullfile() below doesn't double it up
if strcmp(cwd(end), filesep)
    cwd(end) = [];
end

gpnames={'DDx','DDT','DAT','SDx','SDT','SAT'};    %set
pixelSize=0.3642;   % um/px (from .nd2 metadata) -- used for the *_um2 area columns

% Fiber-mask folder/suffix per group.
% [FIX] 舊版只有 4 個元素(來自舊的 4 組資料集),但現在有 12 組,
%       造成 fiberDirs{ggg} 索引超出範圍,且 SDx(第4組)被誤設為 actin_mask_tif。
%       本批 12 組全部都是 a-SMA 染色,所以每一組都用 a-SMA 遮罩。
% fiberDirs    ={'aSMA_mask_tif','aSMA_mask_tif','actin_mask_tif','actin_mask_tif'};
% fiberSuffixes={'_TxRed_aSMAmask.tif','_TxRed_aSMAmask.tif','_TxRed_actinmask.tif','_TxRed_actinmask.tif'};
fiberDirs    =repmat({'aSMA_mask_tif'},1,numel(gpnames));
fiberSuffixes=repmat({'_TxRed_aSMAmask.tif'},1,numel(gpnames));
% [FIX] 防呆:若日後又改了 gpnames 卻忘了更新這兩個陣列,直接報錯提醒
assert(numel(fiberDirs)==numel(gpnames) && numel(fiberSuffixes)==numel(gpnames), ...
    'fiberDirs/fiberSuffixes must have one entry per group.');

%% accumulators for the combined summary (filled in across all groups/images)
AllActinTables={};
AllNucTables={};

%%
for ggg=1:size(gpnames,2)
    %% commonly used folders (built with fullfile so this works on any OS)
    processedDir   = fullfile(cwd, gpnames{ggg}, 'TIFF_files', 'processed');
    binaryDir      = fullfile(processedDir, 'binary_masks');
    actinMaskDir   = fullfile(processedDir, fiberDirs{ggg});
    instanceDir    = fullfile(processedDir, 'instance_masks');
    quantifyDir    = fullfile(processedDir, 'quantification');

    if ~isfolder(processedDir)
        warning('Folder not found - group skipped:\n%s', processedDir);
        continue
    end
    % [FIX] 若該組的 fiber 遮罩資料夾根本不存在(例如第 3 步沒跑到這一組),
    %       一次性提醒並跳過整組,而不是對每個 site 各噴一次警告
    if ~isfolder(actinMaskDir)
        warning('Fiber mask folder not found - group %s skipped (run batch_aSMAFilamentMask for this group first):\n%s', ...
            gpnames{ggg}, actinMaskDir);
        continue
    end
    if ~exist(quantifyDir, 'dir')
        mkdir(quantifyDir);
    end

    %% accumulators for this group's own summary (reset for each group)
    GroupActinTables={};
    GroupNucTables={};

    %% sites are driven by the TxRed SUM files; every other file is built
    % from the same <Well>_S<##> key rather than matched by list position,
    % so a site missing from one step (e.g. failed instance segmentation)
    % is skipped instead of silently pairing images from different sites.
    TxRedFiles=dir(fullfile(processedDir,'*_TxRed_SUM_linear.tif'));
    fprintf('\n=== %s: %d TxRed site(s) ===\n', gpnames{ggg}, numel(TxRedFiles));
    %%
    for ttt=1:size(TxRedFiles,1)
        key=regexprep(TxRedFiles(ttt).name,'_TxRed_SUM_linear\.tif$','');   % e.g. H10_S01
        f.dapi   =fullfile(processedDir,[key '_DAPI_SUM_linear.tif']);
        f.cellBi =fullfile(binaryDir,   [key '_cellBinary.tif']);
        f.nucBi  =fullfile(binaryDir,   [key '_nucleusBinary.tif']);
        f.fiber  =fullfile(actinMaskDir,[key fiberSuffixes{ggg}]);
        f.cellLa =fullfile(instanceDir, [key '_cellLabels.tif']);
        f.nucLa  =fullfile(instanceDir, [key '_nucleusLabels.tif']);
        missing=struct2cell(structfun(@(x) ~isfile(x), f, 'UniformOutput', false));
        if any([missing{:}])
            fn=fieldnames(f);
            warning('%s %s skipped - missing: %s', gpnames{ggg}, key, strjoin(fn([missing{:}]), ', '));
            continue
        end

        %% Real TxRed signals
        disp([gpnames{ggg},' ',key,' TxRed'])
        %
        TxRedImg0=single(imread(fullfile(processedDir,TxRedFiles(ttt).name)));
        TxRedImgBlur=imfilter(TxRedImg0,fspecial('disk',2),'replicate');
        TxRedImgSub=bgsub(TxRedImgBlur,128,0.001);
        %
        TxRedBi=imread(f.cellBi);
        TxRedBgList=TxRedImgSub(TxRedBi==0);
        TxRedBg=GetSeris_bg_bin(TxRedBgList,100);
        TxRedImgReal=TxRedImgSub-TxRedBg;
        %
        figure(100);subplot(1,2,1);imagesc(TxRedImg0);axis image; title('Original');
        figure(100);subplot(1,2,2);imagesc(TxRedImgReal); axis image; title('Real TxRed');

        %% Real DAPI signals
        disp([gpnames{ggg},' ',key,' DAPI'])
        %
        DAPIImg0=single(imread(f.dapi));
        DAPIImgBlur=imfilter(DAPIImg0,fspecial('disk',2),'replicate');
        DAPIImgSub=bgsub(DAPIImgBlur,64,0);
        %
        DAPIBi=imread(f.nucBi);
        DAPIBgList=DAPIImgSub(DAPIBi==0);
        DAPIBg=GetSeris_bg_bin(DAPIBgList,1000);
        DAPIImgReal=DAPIImgSub-DAPIBg;
        %
        figure(200);subplot(1,2,1);imagesc(DAPIImg0);axis image; title('Original');
        figure(200);subplot(1,2,2);imagesc(DAPIImgReal); axis image; title('Real DAPI');

        %% read actin masks and labeled files
        ActinMask=imread(f.fiber)>0;
        TxRedLa=imread(f.cellLa);
        DAPILa=imread(f.nucLa);

        %% quantification of actin objects
        ActinProps=regionprops(ActinMask,'centroid','area','pixelidxlist');
        ActinObjSz=size(ActinProps,1);
        ActinObjMtx=zeros(ActinObjSz,5);
        % 1:x,2:y,3:area,4:CellLabel,5:Intensity
        for aaa=1:ActinObjSz
            ActinObjMtx(aaa,1:2)=ActinProps(aaa).Centroid;
            ActinObjMtx(aaa,3)=ActinProps(aaa).Area;
            TempActinPxIdxList=ActinProps(aaa).PixelIdxList;
            ActinObjMtx(aaa,4)=mode(TxRedLa(TempActinPxIdxList));
            ActinObjMtx(aaa,5)=mean(TxRedImgReal(TempActinPxIdxList));
        end

        %% quantification of nucleated cell objects
        NucLabel0=unique(DAPILa);NucLabel=NucLabel0(NucLabel0>0);
        NucObjSz = numel(NucLabel);
        NucObjMtx = zeros(NucObjSz, 9);
        % 1:Label, 2:CentroidX, 3:CentroidY
        % 4:DAPI-Area, 5:DAPI-Intensity
        % 6:TxRed-Area, 7: TxRed-Intensity
        % 8:Actin-Area, 9: Actin-Intensity
        for bbb = 1:NucObjSz
            % 1:Label, 2:CentroidX, 3:CentroidY
            NucMask = DAPILa == NucLabel(bbb);
            NucProps = regionprops(single(NucMask), 'Centroid', 'Area', 'PixelIdxList');
            NucObjMtx(bbb, 1) = NucLabel(bbb);
            NucObjMtx(bbb, 2:3) = NucProps.Centroid;
            % 4:DAPI-Area, 5:DAPI-Intensity
            NucObjMtx(bbb, 4) = NucProps.Area;
            TempDAPIPxIdxList=NucProps.PixelIdxList;
            NucObjMtx(bbb, 5) = mean(DAPIImgReal(TempDAPIPxIdxList));
            % 6:TxRed-Area, 7: TxRed-Intensity
            CytMask= TxRedLa == NucLabel(bbb);
            CytProps = regionprops(single(CytMask), 'Area','PixelIdxList');
            NucObjMtx(bbb, 6) = CytProps.Area;
            TempTxRedPxIdxList=CytProps.PixelIdxList;
            NucObjMtx(bbb, 7) = mean(TxRedImgReal(TempTxRedPxIdxList));
            % 8:Actin-Area, 9: Actin-Intensity
            CytActinMask= CytMask & ActinMask;
            NucObjMtx(bbb, 8) = sum(CytActinMask, 'all');
            NucObjMtx(bbb, 9) = mean(TxRedImgReal(CytActinMask));

        end

        %% save ActinObjMtx and NucObjMtx as CSV files
        % use the image filename (without extension) as a tag so files
        % from different fields of view don't overwrite each other
        [~, baseName, ~] = fileparts(TxRedFiles(ttt).name);

        ActinTable = array2table(ActinObjMtx, 'VariableNames', ...
            {'CentroidX','CentroidY','Area','CellLabel','Intensity'});
        ActinTable.Area_um2 = ActinTable.Area * pixelSize^2;
        ActinOutFile = fullfile(quantifyDir, [baseName, '_ActinObjMtx.csv']);
        writetable(ActinTable, ActinOutFile);

        NucTable = array2table(NucObjMtx, 'VariableNames', ...
            {'Label','CentroidX','CentroidY','DAPI_Area','DAPI_Intensity', ...
             'TxRed_Area','TxRed_Intensity','Actin_Area','Actin_Intensity'});
        NucTable.DAPI_Area_um2  = NucTable.DAPI_Area  * pixelSize^2;
        NucTable.TxRed_Area_um2 = NucTable.TxRed_Area * pixelSize^2;
        NucTable.Actin_Area_um2 = NucTable.Actin_Area * pixelSize^2;
        NucTable.Actin_Fraction = NucTable.Actin_Area ./ NucTable.TxRed_Area;  % fiber area / cell area
        NucOutFile = fullfile(quantifyDir, [baseName, '_NucObjMtx.csv']);
        writetable(NucTable, NucOutFile);

        %% stash a tagged copy of each table for the group summary and the overall summary
        Group=repmat(gpnames(ggg),size(ActinTable,1),1);
        Image=repmat({baseName},size(ActinTable,1),1);
        TaggedActinTable=[table(Group,Image),ActinTable];
        GroupActinTables{end+1}=TaggedActinTable; %#ok<SAGROW>
        AllActinTables{end+1}=TaggedActinTable; %#ok<SAGROW>

        Group=repmat(gpnames(ggg),size(NucTable,1),1);
        Image=repmat({baseName},size(NucTable,1),1);
        TaggedNucTable=[table(Group,Image),NucTable];
        GroupNucTables{end+1}=TaggedNucTable; %#ok<SAGROW>
        AllNucTables{end+1}=TaggedNucTable; %#ok<SAGROW>

    end

    %% combine this group's per-image tables into one summary file per object type
    if isempty(GroupNucTables)
        warning('No sites quantified in %s - no group summary written.', gpnames{ggg});
        continue
    end
    GroupActinSummary=vertcat(GroupActinTables{:});
    writetable(GroupActinSummary,fullfile(quantifyDir,[gpnames{ggg},'_ActinObjMtx_summary.csv']));

    GroupNucSummary=vertcat(GroupNucTables{:});
    writetable(GroupNucSummary,fullfile(quantifyDir,[gpnames{ggg},'_NucObjMtx_summary.csv']));
end

%% combine all per-image tables into one summary file each, saved at the top level
% [FIX] 若全部組別都沒有量到任何 site,舊版 vertcat 空 cell 後 writetable 會直接報錯
if isempty(AllNucTables)
    error('No sites were quantified in any group - nothing to summarise. Check the "skipped - missing" warnings above.');
end

summaryDir=fullfile(cwd,'quantification_summary');
if ~exist(summaryDir,'dir')
    mkdir(summaryDir);
end

ActinSummary=vertcat(AllActinTables{:});
writetable(ActinSummary,fullfile(summaryDir,'ActinObjMtx_summary.csv'));

NucSummary=vertcat(AllNucTables{:});
writetable(NucSummary,fullfile(summaryDir,'NucObjMtx_summary.csv'));
