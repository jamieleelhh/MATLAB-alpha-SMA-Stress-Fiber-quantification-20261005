%% Step 8 (new, part 3): strong-fiber-positive cells
clear;clc;close all
if ispc
    cwd='E:\浩熏\20260922 HS68(3000) (in SF DMEM + 0.1% BSA) a-SMA Sita(40) TGFb(10) AZD4547 4 days\nd\quantification_summary\';
else
    cwd='/Users/FC/Library/CloudStorage/Dropbox/NTU_Pharm/Projects/2026/LabData/DPP4_20260917/quantification_summary/';
end
gpnames={'DDx','DDT','DAT','SDx','SDT','SAT'};
drugs   ={'D','S'};   drugLbl={'DMSO','Sita'};
conds   ={'Dx','DT','AT'};     condLbl={'-TGFb','DMSO+TGFb','AZD+TGFb'};
cols=[0.55 0.55 0.55; 0.85 0.40 0.10; 0.10 0.40 0.80];

nc=numel(conds);

% ---- strong-fiber threshold ----
refConds  = {'AT'};   % reference groups used to set the threshold
refPct    = 90;       % percentile of reference fiber-object intensities
manualThr = [950];       % e.g. 500 to override the percentile rule

A=readtable([cwd,'ActinObjMtx_summary.csv'],'TextType','char');
N=readtable([cwd,'NucObjMtx_summary.csv'],'TextType','char');
N.Well=cellfun(@(s) s(1:3),N.Image,'UniformOutput',false);
A.Cond=cellfun(@(g) g(end-1:end),A.Group,'UniformOutput',false);

isRef=ismember(A.Cond,refConds) & A.Intensity>0;
thr=prctile(A.Intensity(isRef),refPct);
if ~isempty(manualThr), thr=manualThr; end
fprintf('Strong-fiber threshold = %.1f (a.u.)\n',thr);

%% helpers: does cell (Group|Image|Label) own >=1 fiber object with Intensity >= t ?
keyN=string(N.Group)+"|"+string(N.Image)+"|"+string(N.Label);
keyA=string(A.Group)+"|"+string(A.Image)+"|"+string(A.CellLabel);
strongFlag=@(t) double(ismember(keyN,keyA(A.Intensity>=t)));

[gi,grpI,imgI]=findgroups(N.Group,N.Image);      % per site
[gw,grpW,wellW]=findgroups(N.Group,N.Well);      % per well (pooled cells)
I=table(grpI,imgI,'VariableNames',{'Group','Image'});
W=table(grpW,wellW,'VariableNames',{'Group','Well'});
W.nCells=splitapply(@numel,N.Label,gw);

f=strongFlag(thr);
I.StrongPct=100*splitapply(@mean,f,gi);
W.nStrong=splitapply(@sum,f,gw);
W.StrongPct=100*W.nStrong./W.nCells;
writetable(W,[cwd,'WellLevel_StrongFiberPositive.csv']);   % use for stats (n = wells)
disp(W)

%% Figure 1: superplot at the chosen threshold
figure(100);clf;hold on
for k=1:numel(gpnames)
    g=gpnames{k};
    di=find(strcmp(drugs,g(1:end-2)));  ci=find(strcmp(conds,g(end-1:end)));
    x=(di-1)*(nc+1)+ci;
    yI=I.StrongPct(strcmp(I.Group,g));
    yW=W.StrongPct(strcmp(W.Group,g));
    scatter(x+(rand(size(yI))-0.5)*0.5,yI,14,cols(ci,:),'filled','MarkerFaceAlpha',0.35);
    scatter(x+(rand(size(yW))-0.5)*0.15,yW,80,cols(ci,:),'filled','MarkerEdgeColor','k','LineWidth',1);
    plot(x+[-0.3 0.3],median(yW)*[1 1],'k','LineWidth',2);
end
set(gca,'XTick',(1:numel(drugs))*(nc+1)-(nc+1)/2,'XTickLabel',drugLbl);
xlim([0 numel(drugs)*(nc+1)]); box on
ylabel(sprintf('Cells with >=1 strong fiber (%%)  [thr = %.0f]',thr));
for ci=1:3, scatter(nan,nan,60,cols(ci,:),'filled','DisplayName',condLbl{ci}); end
legend('Location','northeastoutside');
title('small dots = sites, large dots = wells, bar = median of wells');

%% Figure 2: threshold sweep (is the conclusion threshold-dependent?)
thrs=round(linspace(150,1000,18));
PW=zeros(height(W),numel(thrs));
for t=1:numel(thrs)
    ft=strongFlag(thrs(t));
    PW(:,t)=100*splitapply(@sum,ft,gw)./W.nCells;
end
figure(200);clf;tiledlayout(1,numel(drugs),'TileSpacing','compact');
for di=1:numel(drugs)
    nexttile;hold on
    for ci=1:3
        g=[drugs{di},conds{ci}];
        rows=find(strcmp(W.Group,g));
        if isempty(rows), continue; end
        plot(thrs,PW(rows,:)','-','Color',[cols(ci,:) 0.3],'LineWidth',0.8);      % each well
        plot(thrs,mean(PW(rows,:),1),'-','Color',cols(ci,:),'LineWidth',2.5);     % mean of wells
    end
    xline(thr,'k--'); title(drugLbl{di}); xlabel('Strong-fiber threshold (a.u.)');
    if di==1, ylabel('Cells with >=1 strong fiber (%)'); end
    box on
end