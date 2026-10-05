%%
clear;clc;close all
cwd='E:\浩熏\20260922 HS68(3000) (in SF DMEM + 0.1% BSA) a-SMA Sita(40) TGFb(10) AZD4547 4 days\nd\';
%cwd='D:\Documents\NTU_dropbox\Dropbox\NTU_Pharm\Projects\2026\LabData\DPP4_20260917\';
gpnames={'DDx','DDT','DAT','SDx','SDT','SAT'};    %set

showFigures = true;   % set false to run faster (no per-image display)

%%
for ggg=1:size(gpnames,2)
    %%
    procDir = fullfile(cwd, gpnames{ggg}, 'TIFF_files', 'processed');
    if ~isfolder(procDir)
        warning('Folder not found - group skipped:\n%s', procDir);
        continue
    end

    % TxRed and DAPI are looped separately, so a site that is present in one
    % channel but not the other (e.g. skipped during deconvolution) cannot
    % shift the pairing between the two lists.
    TxRedFiles=dir(fullfile(procDir,'*_TxRed_SUM_linear.tif'));
    DAPIFiles=dir(fullfile(procDir,'*_DAPI_SUM_linear.tif'));
    fprintf('=== %s: %d TxRed, %d DAPI ===\n', gpnames{ggg}, numel(TxRedFiles), numel(DAPIFiles));

    %% ---- TxRed ----------------------------------------------------------
    for ttt=1:numel(TxRedFiles)
        disp([gpnames{ggg},' ',num2str(ttt),' TxRed ',TxRedFiles(ttt).name])
        %
        TxRedImg0=single(imread(fullfile(procDir,TxRedFiles(ttt).name)));
        TxRedImgBlur=imfilter(TxRedImg0,fspecial('disk',2),'replicate');
        bg_Blur=GetSeris_bg_bin(TxRedImgBlur,1000);
        TxRedImgSub1=bgsub(TxRedImgBlur,128,0.001);
        TxRedImgSub2=bgsub(TxRedImgBlur,8,0.05);
        TxRedImgSub=(TxRedImgSub1+TxRedImgSub2)/2;
        bg_Sub=GetSeris_bg_bin(TxRedImgSub,1000);
        TxRedImgFinal=uint16(TxRedImgSub+(bg_Blur-bg_Sub));

        %%
        if showFigures
            figure(100);subplot(1,2,1);imagesc(TxRedImg0);axis image; title(TxRedFiles(ttt).name,'Interpreter','none');
            figure(100);subplot(1,2,2);imagesc(TxRedImgFinal); axis image; title('Background Subtracted Image');
            drawnow
        end

        %%
        imwrite(TxRedImgFinal,fullfile(procDir,[TxRedFiles(ttt).name(1:end-4),'_BgSubtracted.tif']), ...
            'compression','none')
    end

    %% ---- DAPI -----------------------------------------------------------
    for ttt=1:numel(DAPIFiles)
        disp([gpnames{ggg},' ',num2str(ttt),' DAPI ',DAPIFiles(ttt).name])
        %
        DAPIImg0=single(imread(fullfile(procDir,DAPIFiles(ttt).name)));
        DAPIImgBlur=imfilter(DAPIImg0,fspecial('disk',2),'replicate');
        bg_Blur=GetSeris_bg_bin(DAPIImgBlur,1000);
        DAPIImgSub1=bgsub(DAPIImgBlur,64,0.01);
        DAPIImgSub2=bgsub(DAPIImgBlur,16,0.05);
        DAPIImgSub=(DAPIImgSub1+DAPIImgSub2)/2;
        bg_Sub=GetSeris_bg_bin(DAPIImgSub,1000);
        DAPIImgFinal=uint16(DAPIImgSub+(bg_Blur-bg_Sub));

        %%
        if showFigures
            figure(200);subplot(1,2,1);imagesc(DAPIImg0);axis image; title(DAPIFiles(ttt).name,'Interpreter','none');
            figure(200);subplot(1,2,2);imagesc(DAPIImgFinal); axis image; title('Background Subtracted Image');
            drawnow
        end

        %%
        imwrite(DAPIImgFinal,fullfile(procDir,[DAPIFiles(ttt).name(1:end-4),'_BgSubtracted.tif']), ...
            'compression','none')
    end
end
