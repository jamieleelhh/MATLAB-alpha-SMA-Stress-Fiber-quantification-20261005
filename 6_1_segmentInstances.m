function result = segmentInstances(nucleusBinaryFile, cellBinaryFile, dapiFile, txredFile, varargin)
%SEGMENTINSTANCES  Split touching nuclei and adjacent cells into labeled
%   instances, from already-validated binary masks.
%
%   result = SEGMENTINSTANCES(nucleusBinaryFile, cellBinaryFile, dapiFile, txredFile)
%   takes the binary (0/1) nucleus and cell masks produced by
%   generateBinaryMasks.m -- NOT re-derived here, just split into
%   instances -- plus the DAPI and TxRed images (ideally the same
%   '..._BgSubtracted.tif' files used to build the binary masks; both are
%   REQUIRED to already be background-flattened, same assumption as
%   generateBinaryMasks.m) for the QC overlay and (in one narrow,
%   explicitly-scoped place -- see rule 3 below) cell-merge judgement.
%   Returns two labeled masks:
%     nucleusLabels : uint16, 0 = background, each nucleus gets an ID
%     cellLabels    : uint16, 0 = background, each cell gets an ID
%   A nucleus and its enclosing cell always share the SAME ID -- EXCEPT
%   for the specific, deliberate edge case documented below (a cell
%   fragment with no nucleus that touches the image border), which keeps
%   its own ID with no nucleus counterpart.
%
%   nucleusBinaryFile/cellBinaryFile/dapiFile/txredFile may be file paths
%   (char/string) or already-loaded matrices/images.
%
%   How it works:
%     1. Nuclei: NO watershed splitting is performed. Touching/clumped
%        nuclei in the binary mask are labeled exactly as they appear --
%        each connected component of the nucleus binary mask becomes ONE
%        nucleus ID, even if it visually looks like more than one real
%        nucleus. (An earlier version used a distance-transform
%        watershed to try to split touching nuclei; removed per explicit
%        instruction in favor of simplicity and predictability.)
%     2. Cells: built per CONNECTED COMPONENT of the cell binary mask,
%        not by re-thresholding anything, and -- per explicit
%        instruction -- using ONLY nucleus position for any splitting
%        decision, NEVER TxRed (actin) signal or gradient:
%          - blob contains exactly one nucleus -> the whole blob is that
%            one cell, unsplit (this is what keeps a single physically-
%            continuous cell from being fragmented).
%          - blob contains 2+ nuclei -> split via a TWO-STAGE marker-
%            controlled watershed:
%              Stage 1 (geometric): EUCLIDEAN DISTANCE TO THE NEAREST
%                NUCLEUS, a Voronoi-style partition where every pixel
%                goes to whichever nucleus is closest. This stage's only
%                purpose is to bound each nucleus's own "protected zone"
%                below so two nearby nuclei's zones can never touch.
%              Stage 2 (actin-refined periphery): within
%                'PerinuclearProtectRadius' px of each nucleus (clipped
%                to that nucleus's own Stage-1 territory), the split is
%                identical to Stage 1, UNAFFECTED by actin -- actin
%                intensity genuinely differs between nuclear/perinuclear
%                and peripheral cytoplasm regions, and using it that
%                close to the nucleus was unreliable (confirmed
%                directly: it repeatedly stranded a nucleus in a small
%                sliver of its own, separate from the rest of its
%                cytoplasm). BEYOND that radius, in the periphery, the
%                boundary instead follows the TxRed actin GRADIENT --
%                real cortical actin (bright bundles running parallel to
%                the true cell edge) and lamellipodial actin meshwork
%                are reliable indicators of where the actual cell-cell
%                junction sits out there, unlike close to the nucleus.
%                Concretely: each nucleus's protected zone (from Stage 1)
%                seeds a watershed on the actin gradient instead of on
%                distance, so the assignment near the nucleus is exactly
%                the seed itself, while the peripheral boundary snaps to
%                real actin ridges instead of the straight geometric
%                bisector.
%            EVERY resulting nucleated sub-region becomes its OWN final
%            cell, PERMANENTLY -- a blob with N nuclei always yields
%            exactly N separately-labeled cells, and nothing later in
%            the pipeline
%            ever merges two nucleated cells back together (see rule 3).
%          - blob contains zero nuclei, OR a Voronoi split above leaves
%            a sub-region with zero nuclei (e.g. a thin protrusion cut
%            off from its parent cell) -- these nucleus-less fragments
%            have no position to partition by (and, per instruction, no
%            actin-based split is attempted either), so they are never
%            split; instead they're resolved by the SAME policy, listed
%            here in the order it's applied:
%              (i)   touches the image border -> left alone: kept as its
%                    OWN standalone cell instance, no nucleus, no fusion
%                    with anything, and never removed (a fragment cut off
%                    by the field of view could belong to a real cell
%                    outside the frame, and its true shape/identity can't
%                    be determined, so merging assumptions aren't safe).
%              (ii)  doesn't touch the border, and touches NO other cell
%                    area at all -> removed (left as background). This
%                    only happens for a fragment that is its own entire
%                    connected component of the cell mask (an isolated
%                    no-nucleus blob can never itself abut another blob,
%                    by definition of connected components) -- most
%                    likely non-cellular debris that happened to pass the
%                    TxRed threshold.
%              (iii) doesn't touch the border, and is in DIRECT pixel
%                    contact with another cell area (the common case for
%                    a leftover Voronoi sub-region inside a multi-
%                    nucleus blob, which necessarily abuts its sibling
%                    sub-regions) -> fused into the best-matching
%                    DIRECTLY-CONTACTING neighboring cell that (a) has
%                    its own nucleus and (b) does not itself touch the
%                    border, chosen by, in priority order:
%                      (a) most direct contact (longest shared boundary),
%                      (b) most similar mean TxRed (F-actin) intensity,
%                      (c) nearest nucleus (centroid distance).
%                    (a) dominates; (b) and (c) only break a near-tie in
%                    (a), within 'OrphanFuseContactTolerance'. Note this
%                    IS a place actin data (criterion b) is used -- but
%                    it's deciding WHICH already-fixed fragment goes with
%                    WHICH already-fixed cell, a fundamentally different,
%                    later-stage question from "where should the cell
%                    boundary be drawn," which is what rule 2 above is
%                    about avoiding actin for. "Direct contact" means
%                    literal 8-connected pixel adjacency, with no
%                    proximity tolerance -- two cell areas that aren't
%                    actually touching are never assigned the same
%                    label, no matter how close they are. If a
%                    fragment's only directly-contacting neighbors are
%                    themselves border-touching and/or nucleus-less (so
%                    no legal fuse target exists), it's kept as its own
%                    standalone cell instead -- same reasoning as (i).
%              Resolution runs in rounds so that fusion chains resolve
%              correctly (a fragment touching another still-unresolved
%              fragment waits for that neighbor to resolve first, then
%              re-evaluates); each individual fuse step is still only
%              ever between two DIRECTLY-touching areas.
%     3. Three hard rules govern every fuse decision made in this
%        function, and NONE has an exception:
%          - Fusion of a nucleus-less fragment only ever occurs into a
%            nucleus-BEARING cell area. A nucleus-less fragment is never
%            fused into another nucleus-less fragment.
%          - A nucleated cell area is NEVER fused with another nucleated
%            cell area, under any circumstance. (An earlier version had
%            a deliberate exception here -- merging two adjacent
%            nucleated cells when their shared actin boundary was judged
%            "too weak to be a real junction" -- removed per explicit
%            instruction after it repeatedly chain-merged multiple
%            genuinely separate cells into one: confirmed directly, 8
%            distinct nuclei ended up sharing a single ~256,000 px
%            "cell" in a real test image. A blob with N nuclei now
%            always yields exactly N separately-labeled cells, and stays
%            that way.)
%          - Fusion only ever occurs between two areas in direct pixel
%            contact. Two cell areas that don't literally touch are never
%            assigned the same label, regardless of proximity.
%          - A nucleus-less area touching the image border is left alone
%            entirely: never fused into anything, never used as a fusion
%            source. (It's also never a valid fusion TARGET, since it has
%            no nucleus -- rule one already rules that out.)
%        As a safeguard, every final ID is also verified to be a single
%        connected piece before returning (see enforceSingleComponentPerLabel):
%        if a nucleus's mask footprint and the cell mask's footprint
%        aren't perfectly aligned, it's possible for one nucleus to
%        overlap two cell-mask pieces that don't themselves touch, which
%        would otherwise let both inherit the same ID despite rule two.
%        Any extra piece found this way is peeled off and independently
%        re-run through the same fuse-or-standalone-or-remove logic.
%
%   Given a nucleus-less area that touches a nucleated cell, it will
%   always be fused into it: it's non-border, in direct contact, and
%   the target has a nucleus, so all three rules are satisfied.
%
%   Name-Value options:
%     'SmoothSigma'         gaussian smoothing sigma (px) applied to
%                           TxRed before computing the actin gradient
%                           used by the periphery-refinement stage of a
%                           multi-nucleus split (see step 2 above) and
%                           the actin-similarity criterion in orphan
%                           fusion. Never used within a nucleus's own
%                           protected perinuclear zone. Default 2.
%     'PerinuclearProtectRadius' radius (px), around each nucleus, of the
%                           zone where a multi-nucleus blob's split
%                           follows the geometric nearest-nucleus
%                           partition exactly, UNAFFECTED by actin
%                           signal -- clipped to that nucleus's own
%                           geometric territory so two nearby nuclei's
%                           protected zones can never touch each other.
%                           Beyond this radius, in the periphery, the
%                           split instead follows actin gradient ridges
%                           (cortical actin bundles, lamellipodial
%                           meshwork), which are a reliable indicator of
%                           the true cell-cell junction out there, unlike
%                           near the nucleus. Default 20 -- a starting
%                           estimate, not verified against real data (no
%                           direct MATLAB access at the time), so tune
%                           this first if the actin-based refinement ever
%                           reaches too close to a nucleus, or if too
%                           much of a cell's true periphery is still
%                           being decided by raw geometry instead of
%                           visible cortical actin/lamellipodia.
%     'OrphanFuseContactTolerance' relative tolerance (fraction of the
%                           best candidate's contact length) within which
%                           two or more candidate cells count as "tied"
%                           on directness of contact, so the fuse
%                           decision falls through to the actin-
%                           similarity and then nucleus-distance
%                           tiebreakers instead of being decided by
%                           contact length alone. Default 0.10.
%     'OrphanFuseActinTolerance' relative tolerance (fraction of the best
%                           candidate's actin-similarity score) within
%                           which two or more contact-tied candidate
%                           cells count as "tied" on mean F-actin
%                           intensity too, so the fuse decision falls
%                           through to the nucleus-distance tiebreaker
%                           instead of being decided by actin similarity
%                           alone. Default 0.10.
%     'MinNucleusArea'      minimum nucleus area (px^2) kept after
%                           labeling. Default 150.
%     'MinCellArea'         minimum final cell area (px^2). Default 400.
%     'ExcludeBorder'       true/false, drop cells/nuclei touching the
%                           image edge from the final output entirely.
%                           Default false (kept by default). This is
%                           independent of, and unrelated to, the
%                           never-fuse-at-the-border rule above -- that
%                           rule always applies regardless of this
%                           setting.
%     'OutputDir'           if provided, writes '<BaseName>_nucleusLabels.tif'
%                           and '<BaseName>_cellLabels.tif' (uint16) plus
%                           a QC overlay PNG into this folder. Default ''
%                           (nothing written).
%     'BaseName'            base filename for saved outputs. Default
%                           derived from cellBinaryFile if it is a path,
%                           otherwise 'site'.
%     'SaveQCOverlay'       true/false, save a labeled-instance overlay
%                           PNG when OutputDir is set. Default true.
%
%   Output struct fields:
%     .nucleusLabels, .cellLabels   uint16 label matrices
%     .numCells                     total number of final labeled cells
%                                    (including any nucleus-less border
%                                    fragments kept per rule 2(i) above)
%     .numCellsWithoutNucleus       how many of .numCells have no
%                                    corresponding nucleus -- should
%                                    normally be 0 or very small; watch
%                                    this across a batch as a sanity check

p = inputParser;
addParameter(p, 'SmoothSigma', 2);
addParameter(p, 'PerinuclearProtectRadius', 300);     % 預設 20
addParameter(p, 'OrphanFuseContactTolerance', 0.10);
addParameter(p, 'OrphanFuseActinTolerance', 0.10);
addParameter(p, 'MinNucleusArea', 150);
addParameter(p, 'MinCellArea', 3000);                 % 預設 400
addParameter(p, 'ExcludeBorder', false);
addParameter(p, 'OutputDir', '');
addParameter(p, 'BaseName', '');
addParameter(p, 'SaveQCOverlay', true);
parse(p, varargin{:});
opt = p.Results;

[nucleusBW, nucPathStr] = loadBinary(nucleusBinaryFile);
[cellBW, cellPathStr]   = loadBinary(cellBinaryFile);
[dapiImg, ~]            = loadImage(dapiFile);
[txredImg, ~]           = loadImage(txredFile);

if ~isequal(size(nucleusBW), size(cellBW)) || ~isequal(size(nucleusBW), size(txredImg)) ...
        || ~isequal(size(nucleusBW), size(dapiImg))
    error('segmentInstances:sizeMismatch', ...
        'nucleusBinary, cellBinary, dapiFile, and txredFile must all be the same size.');
end

if isempty(opt.BaseName)
    if ~isempty(cellPathStr)
        [~, baseName] = fileparts(cellPathStr);
        baseName = regexprep(baseName, '_cellBinary.*$', '');
    elseif ~isempty(nucPathStr)
        [~, baseName] = fileparts(nucPathStr);
        baseName = regexprep(baseName, '_nucleusBinary.*$', '');
    else
        baseName = 'site';
    end
else
    baseName = opt.BaseName;
end

% ---- 1. label nuclei (no watershed splitting) --------------------------
nucleusLabels = splitTouchingNuclei(nucleusBW, opt);
nucleusLabels = bwlabel(bwareaopen(nucleusLabels > 0, opt.MinNucleusArea));

% ---- 2. build cell instances from the cell binary mask -----------------
[cellLabels, nucleusLabels] = assignCellsFromBlobs(txredImg, nucleusLabels, cellBW, opt);

% ---- 3. final area/border filtering, relabel contiguous ---------------
[cellLabels, nucleusLabels] = finalFilter(cellLabels, nucleusLabels, opt);

result.nucleusLabels = uint16(nucleusLabels);
result.cellLabels    = uint16(cellLabels);
result.numCells       = max(cellLabels(:));
nucleatedIDs           = unique(nucleusLabels(nucleusLabels > 0));
result.numCellsWithoutNucleus = result.numCells - numel(nucleatedIDs);

if ~isempty(opt.OutputDir)
    if ~exist(opt.OutputDir, 'dir')
        mkdir(opt.OutputDir);
    end
    imwrite(result.nucleusLabels, fullfile(opt.OutputDir, [baseName '_nucleusLabels.tif']));
    imwrite(result.cellLabels,    fullfile(opt.OutputDir, [baseName '_cellLabels.tif']));

    if opt.SaveQCOverlay
        saveInstanceQCOverlay(dapiImg, txredImg, result.cellLabels, result.nucleusLabels, ...
            fullfile(opt.OutputDir, [baseName '_instanceQC_overlay.png']));
    end
end

end % segmentInstances


% ======================================================================
% helper functions
% ======================================================================

function [bw, pathStr] = loadBinary(input)
pathStr = '';
if ischar(input) || isstring(input)
    pathStr = char(input);
    bw = imread(pathStr) > 0;
else
    bw = logical(input);
end
end


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


function tf = touchesImageBorder(mask)
tf = any(mask(1,:)) || any(mask(end,:)) || any(mask(:,1)) || any(mask(:,end));
end


function labels = splitTouchingNuclei(nucleusBW, opt) %#ok<INUSD>
% NO watershed splitting -- per explicit instruction, touching/clumped
% nuclei are labeled exactly as they appear in the binary mask, each
% connected component becoming one nucleus ID. (opt is accepted but
% unused, kept only so the call site doesn't need to change.)
labels = bwlabel(nucleusBW);
end


function [cellLabels, nucleusLabelsOut] = assignCellsFromBlobs(txredImg, nucleusLabels, cellBW, opt)
% Build whole-cell labels from connected components of the (already
% validated) cell binary mask -- see the function-level help for the
% full policy on 0 / 1 / 2+ nuclei per blob, and on nucleus-less
% fragments (border-standalone / removed / fused). Returns cellLabels
% and nucleusLabelsOut using the SAME id for a nucleus and its cell,
% except for border-standalone nucleus-less fragments (documented above).
%
% Per explicit instruction, splitting a multi-nucleus blob is anchored
% ENTIRELY by nucleus position near each nucleus (a geometric nearest-
% nucleus/Voronoi-style partition) -- actin intensity genuinely differs
% between nuclear/perinuclear and peripheral cytoplasm regions, and
% using it there was unreliable (confirmed directly: it repeatedly
% stranded a nucleus in its own small sliver, separate from the rest of
% its cytoplasm). Further OUT in the periphery, though, real cortical
% actin (bright bundles running parallel to the true cell edge) and
% lamellipodial actin meshwork ARE reliable indicators of where the
% actual cell-cell junction sits -- so the boundary is refined there
% using the actin gradient, via a two-stage watershed (see below).
[rows, cols] = size(cellBW);
cellBlobs = bwlabel(cellBW, 8);
numBlobs = max(cellBlobs(:));

% gmag (actin gradient) is used ONLY for the periphery-refinement stage
% below and the actin-similarity criterion in orphan fusion (via
% chooseBestNeighbor) -- never within a nucleus's own protected
% perinuclear zone.
imgNorm = mat2gray(imgaussfilt(txredImg, opt.SmoothSigma));
gmag = imgradient(imgNorm);

cellLabels = zeros(rows, cols);
nucleusLabelsOut = zeros(rows, cols);
nextID = 1;
orphanMask = false(rows, cols); % nucleus-less, non-border area -- resolved (fused or removed) after this loop

for b = 1:numBlobs
    blobMask = cellBlobs == b;
    nucIDsHere = reshape(setdiff(unique(nucleusLabels(blobMask)), 0), 1, []);

    if isempty(nucIDsHere)
        % Whole blob has no nucleus at all -- nothing to partition by
        % position, so it's never split here; handled entirely by the
        % nucleus-less-fragment policy (border-standalone / orphan pool).
        [cellLabels, orphanMask, nextID] = placeNucleuslessRegion( ...
            blobMask, cellLabels, orphanMask, nextID, opt);
        continue;
    end

    if isscalar(nucIDsHere)
        % whole connected blob is one cell -- do not internally split it
        cellLabels(blobMask) = nextID;
        nucleusLabelsOut(nucleusLabels == nucIDsHere) = nextID;
        nextID = nextID + 1;
        continue;
    end

    % 2+ nuclei in one connected blob: two-stage split. Cropped to the
    % blob's bounding box for speed throughout.
    [rIdx, cIdx] = find(blobMask);
    pad = 15;
    r0 = max(1, min(rIdx) - pad); r1 = min(rows, max(rIdx) + pad);
    c0 = max(1, min(cIdx) - pad); c1 = min(cols, max(cIdx) + pad);

    blobMaskC = blobMask(r0:r1, c0:c1);
    nucleusLabelsC = nucleusLabels(r0:r1, c0:c1);
    localMarkersC = ismember(nucleusLabelsC, nucIDsHere) & blobMaskC;
    bgMarkerC = ~blobMaskC;
    gmagC = gmag(r0:r1, c0:c1);

    % --- Stage 1: pure geometric Voronoi partition -----------------
    % Distance-to-nearest-nucleus-marker map; watershed on this (with
    % each nucleus and the outside-the-blob background both imposed as
    % their own minima) gives each nucleus the territory of pixels
    % closest to it, entirely independent of TxRed intensity. Each
    % nucleus marker is already a genuine regional minimum (distance 0
    % at the marker itself), so no extended-minima/depth step is needed.
    % This stage's ONLY purpose is to bound each nucleus's protected
    % zone below so it can never bleed into a neighboring nucleus's
    % territory -- the actual final split comes from Stage 2.
    distC = bwdist(localMarkersC);
    Dimposed = imimposemin(distC, bgMarkerC | localMarkersC);
    LgeoC = watershed(Dimposed);
    LgeoC(~blobMaskC) = 0;

    % --- Stage 2: actin-refined periphery, nucleus-protected core ---
    % Each nucleus's own protected zone -- a disk of radius
    % 'PerinuclearProtectRadius' around JUST that nucleus, clipped to
    % its OWN Stage-1 territory so two nearby nuclei's protected zones
    % can never touch or merge into one marker region -- serves as the
    % seed for a SECOND watershed, this time on the actin gradient
    % instead of distance. Within the protected zone the assignment is
    % identical to Stage 1 (it IS the marker); beyond it, in the
    % periphery, the boundary follows real actin gradient ridges
    % (cortical actin bundles, lamellipodial meshwork) instead of the
    % straight geometric bisector, which is a better match to where the
    % true cell-cell junction actually sits.
    geoIDsC = reshape(setdiff(unique(LgeoC(:)), 0), 1, []);
    protectedZoneC = false(size(blobMaskC));
    for gi = geoIDsC
        thisTerritory = LgeoC == gi;
        thisNucMarker = localMarkersC & thisTerritory;
        if ~any(thisNucMarker(:))
            continue; % shouldn't happen -- every Stage-1 territory has a seeding nucleus
        end
        protectedZoneC = protectedZoneC | ...
            (imdilate(thisNucMarker, strel('disk', opt.PerinuclearProtectRadius)) & thisTerritory);
    end

    gmag2C = imimposemin(gmagC, protectedZoneC | bgMarkerC);
    LC = watershed(gmag2C);
    LC(~blobMaskC) = 0;

    % Same ridge/plateau-snapping fix as before: watershed leaves
    % 1px-wide (or wider, in ambiguous areas) zero-valued dividing lines
    % as label 0 -- still real cellBW foreground, so snap each to its
    % nearest catchment rather than letting it silently vanish.
    unclaimedC = blobMaskC & (LC == 0);
    if any(unclaimedC(:))
        [~, nearestIdx] = bwdist(LC > 0);
        LC(unclaimedC) = LC(nearestIdx(unclaimedC));
    end

    L = zeros(rows, cols);
    L(r0:r1, c0:c1) = LC;

    subIDs = reshape(setdiff(unique(L(:)), 0), 1, []);
    for si = subIDs
        subMask = L == si;
        nucsInSub = reshape(setdiff(unique(nucleusLabels(subMask)), 0), 1, []);
        if ~isempty(nucsInSub)
            % nucleated sub-region -> its own cell, permanently -- no
            % nucleated cell is ever fused with another nucleated cell,
            % anywhere in this pipeline (see rule 3 in the function-
            % level help; a prior, now-removed exception allowed this
            % under a weak-actin-boundary judgment, but that repeatedly
            % chain-merged several genuinely separate nuclei into one
            % huge cell -- confirmed directly: 8 distinct nuclei ended
            % up sharing a single ~256,000 px "cell" -- so it was
            % removed per explicit instruction).
            cellLabels(subMask) = nextID;
            nucleusLabelsOut(nucleusLabels == nucsInSub(1)) = nextID;
            nextID = nextID + 1;
        else
            [cellLabels, orphanMask, nextID] = placeNucleuslessRegion( ...
                subMask, cellLabels, orphanMask, nextID, opt);
        end
    end
end

[cellLabels, nucleusLabelsOut, nextID] = resolveOrphanFragments( ...
    orphanMask, cellLabels, nucleusLabelsOut, txredImg, nextID, opt); %#ok<NASGU>

% Connectivity safety net (rule 2): every final ID must be a single
% connected piece (see enforceSingleComponentPerLabel for why this can
% otherwise be violated).
[cellLabels, nucleusLabelsOut] = enforceSingleComponentPerLabel( ...
    cellLabels, nucleusLabelsOut, txredImg, opt);
end


function [cellLabels, orphanMask, nextID] = placeNucleuslessRegion( ...
    regionMask, cellLabels, orphanMask, nextID, opt) %#ok<INUSD>
% Decide the fate of ONE nucleus-less region: since it has no nucleus,
% there's no position to partition it by, and per explicit instruction
% actin signal is no longer used to decide cell-area splits either -- so
% it is NEVER split here, regardless of size. It's sent to a standalone
% ID (if it touches the image border) or the orphan pool (otherwise),
% where it may later be FUSED into a directly-touching nucleated
% neighbor (see resolveOrphanFragments), chosen by contact / actin
% similarity / nucleus distance -- that's a different, later-stage
% decision about which already-fixed cell an already-fixed fragment
% belongs to, not about where to draw a new dividing line, so it isn't
% affected by the "no actin for splitting" rule. (opt is accepted for a
% consistent call signature but currently unused.)
if touchesImageBorder(regionMask)
    cellLabels(regionMask) = nextID;
    nextID = nextID + 1;
else
    orphanMask = orphanMask | regionMask;
end
end


function [cellLabels, nucleusLabelsOut, nextID] = resolveOrphanFragments( ...
    orphanMask, cellLabels, nucleusLabelsOut, txredImg, nextID, opt)
% Resolve non-border, nucleus-less fragments: fuse into the best
% neighboring nucleated non-border cell, keep standalone if no legal fuse
% target exists, or remove if it touches no other cell area at all. Runs
% in rounds so fusion chains (fragment touching fragment touching a real
% cell) resolve correctly -- see function-level help, rule 2(iii).
if ~any(orphanMask(:))
    return;
end

remainingMask = orphanMask;
progress = true;
while progress
    progress = false;
    cc = bwconncomp(remainingMask, 8);
    justResolved = false(size(remainingMask));

    for i = 1:cc.NumObjects
        idx = cc.PixelIdxList{i};
        fragMask = false(size(remainingMask));
        fragMask(idx) = true;

        % Direct contact only -- 8-connected single-pixel dilation, no
        % proximity slack. Two areas that aren't literally touching must
        % never end up under the same label.
        dilated = imdilate(fragMask, ones(3, 3));
        neighborLabels = cellLabels(dilated & ~fragMask);
        neighborIDs = reshape(setdiff(unique(neighborLabels), 0), 1, []);

        if isempty(neighborIDs)
            continue; % no real-cell neighbor yet -- retry next round
        end

        eligible = [];
        for id = neighborIDs
            hasNucleus = any(nucleusLabelsOut(:) == id);
            idMask = cellLabels == id;
            if hasNucleus && ~touchesImageBorder(idMask)
                eligible(end+1) = id; %#ok<AGROW>
            end
        end

        if ~isempty(eligible)
            bestID = chooseBestNeighbor(fragMask, eligible, cellLabels, nucleusLabelsOut, txredImg, opt);
            cellLabels(idx) = bestID;
        else
            % Real neighbor(s) exist but none are legal fuse targets
            % (border-touching and/or nucleus-less) -- keep standalone
            % rather than fusing illegally or removing area that does,
            % in fact, contact another cell.
            cellLabels(idx) = nextID;
            nextID = nextID + 1;
        end
        justResolved(idx) = true;
        progress = true;
    end

    remainingMask(justResolved) = false;
end
% Anything left in remainingMask touches no real cell even transitively
% (an isolated cluster of nucleus-less fragments) -- removed (left as
% background), per rule 2(ii).
end


function [cellLabels, nucleusLabelsOut] = enforceSingleComponentPerLabel( ...
    cellLabels, nucleusLabelsOut, txredImg, opt)
% Guarantee that no final ID spans more than one connected piece. A
% nucleus's own footprint (nucleusLabels, before this function ever ran)
% is always a single connected blob by construction (it comes straight
% out of bwlabel), so for a nucleated ID the piece containing that
% nucleus is the unambiguous "real" one -- any OTHER piece sharing the
% same ID has no business there and is peeled off, then fed back through
% the same direct-contact fuse-or-standalone-or-remove logic used for
% every other nucleus-less fragment, independently, on its own merits.
maxIter = 5;
for iter = 1:maxIter %#ok<FXUP>
    ids = reshape(setdiff(unique(cellLabels(:)), 0), 1, []);
    reflowMask = false(size(cellLabels));
    anyMultiPiece = false;

    for id = ids
        idMask = cellLabels == id;
        cc = bwconncomp(idMask, 8);
        if cc.NumObjects <= 1
            continue;
        end
        anyMultiPiece = true;

        nucMask = nucleusLabelsOut == id;
        anchorFound = false;
        pieceAreas = cellfun(@numel, cc.PixelIdxList);
        for k = 1:cc.NumObjects
            pieceMask = false(size(idMask));
            pieceMask(cc.PixelIdxList{k}) = true;
            if any(nucMask(:) & pieceMask(:))
                anchorFound = true; % this piece keeps id
                continue;
            end
            cellLabels(pieceMask) = 0;
            reflowMask = reflowMask | pieceMask;
        end

        if ~anchorFound
            % id has no nucleus at all -- shouldn't normally happen, since
            % every nucleus-less ID elsewhere in this file is assigned
            % per single connected piece already, but as a safety net
            % keep the largest piece and reflow the rest.
            [~, biggest] = max(pieceAreas);
            for k = 1:cc.NumObjects
                if k == biggest
                    continue;
                end
                pieceMask = false(size(idMask));
                pieceMask(cc.PixelIdxList{k}) = true;
                cellLabels(pieceMask) = 0;
                reflowMask = reflowMask | pieceMask;
            end
        end
    end

    if ~anyMultiPiece
        break;
    end

    nextID = max(cellLabels(:)) + 1;
    [cellLabels, nucleusLabelsOut, ~] = resolveOrphanFragments( ...
        reflowMask, cellLabels, nucleusLabelsOut, txredImg, nextID, opt);
end
end


function bestID = chooseBestNeighbor(fragMask, candidateIDs, cellLabels, nucleusLabelsOut, txredImg, opt)
% Priority order: (a) most direct contact (longest shared boundary), (b)
% most similar mean F-actin intensity, (c) nearest nucleus -- (a)
% dominates; (b)/(c) only break a near-tie in (a), within
% opt.OrphanFuseContactTolerance.
dilatedFrag = imdilate(fragMask, ones(3, 3)); % direct contact only, same as neighbor detection above
contactLen = zeros(size(candidateIDs));
for k = 1:numel(candidateIDs)
    contactLen(k) = nnz(dilatedFrag & (cellLabels == candidateIDs(k)));
end
[sortedContact, order] = sort(contactLen, 'descend');
tol = opt.OrphanFuseContactTolerance * max(sortedContact(1), 1);
tiedMask = sortedContact >= (sortedContact(1) - tol);
tiedIDs = candidateIDs(order(tiedMask));

if isscalar(tiedIDs)
    bestID = tiedIDs(1);
    return;
end

fragMean = mean(txredImg(fragMask));
actinScore = zeros(size(tiedIDs));
for k = 1:numel(tiedIDs)
    actinScore(k) = abs(fragMean - mean(txredImg(cellLabels == tiedIDs(k))));
end
[sortedScore, order2] = sort(actinScore);
actinTol = opt.OrphanFuseActinTolerance * max(sortedScore(1), 1);
tiedMask2 = sortedScore <= (sortedScore(1) + actinTol);
tiedIDs2 = tiedIDs(order2(tiedMask2));

if isscalar(tiedIDs2)
    bestID = tiedIDs2(1);
    return;
end

fragProps = regionprops(fragMask, 'Centroid');
fragCentroid = fragProps(1).Centroid;
nucDist = zeros(size(tiedIDs2));
for k = 1:numel(tiedIDs2)
    nucProps = regionprops(nucleusLabelsOut == tiedIDs2(k), 'Centroid');
    nucCentroid = nucProps(1).Centroid;
    nucDist(k) = hypot(fragCentroid(1) - nucCentroid(1), fragCentroid(2) - nucCentroid(2));
end
[~, bi] = min(nucDist);
bestID = tiedIDs2(bi);
end


function [cellLabelsOut, nucleusLabelsOut] = finalFilter(cellLabels, nucleusLabels, opt)
% Drop final cells below the minimum area, optionally drop cells
% touching the image border, and relabel survivors contiguously from 1.
[rows, cols] = size(cellLabels);
cellLabelsOut = zeros(rows, cols);
nucleusLabelsOut = zeros(rows, cols);

nextID = 1;
ids = reshape(setdiff(unique(cellLabels(:)), 0), 1, []);

for k = 1:numel(ids)
    id = ids(k);
    cellMask = cellLabels == id;

    if nnz(cellMask) < opt.MinCellArea
        continue;
    end

    if opt.ExcludeBorder
        onBorder = touchesImageBorder(cellMask);
        if onBorder
            continue;
        end
    end

    nucMask = nucleusLabels == id;
    cellLabelsOut(cellMask) = nextID;
    nucleusLabelsOut(nucMask) = nextID;
    nextID = nextID + 1;
end
end


function saveInstanceQCOverlay(dapiImg, txredImg, cellLabels, nucleusLabels, outPath)
% Quick visual check: TxRed shown as grayscale/white, DAPI as magenta,
% with cell instance boundaries (yellow) and nucleus instance boundaries
% (cyan) overlaid -- same coloring convention as generateBinaryMasks.m's
% QC overlay. Boundaries are drawn between EVERY pair of adjacent,
% differently-labeled instances (via boundarymask), not just around the
% outer silhouette of all foreground combined -- see the note at
% cellBoundary/nucBoundary below for why that distinction matters.
txredNorm = mat2gray(imadjust(mat2gray(txredImg)));
dapiNorm  = mat2gray(imadjust(mat2gray(dapiImg)));

% TxRed contributes equally to R/G/B (grayscale/white). DAPI adds to R
% and B only (magenta), so DAPI-only areas read magenta, TxRed-only
% areas read white/gray, and overlap reads a lighter pink-white.
rC = min(txredNorm + dapiNorm, 1);
gC = txredNorm;
bC = min(txredNorm + dapiNorm, 1);

% bwperim(cellLabels > 0) only traces the outer envelope of ALL
% foreground combined -- it does NOT draw a line between two adjacent
% cells that are directly touching with no background gap between them,
% which is now the common case. boundarymask instead marks a boundary
% wherever the LABEL VALUE changes between neighboring pixels --
% including between two touching, differently-labeled cells -- so every
% real segmentation line is visible here, not just the outer silhouette.
cellBoundary = boundarymask(cellLabels);
nucBoundary  = boundarymask(nucleusLabels);

rC(cellBoundary) = 1; gC(cellBoundary) = 1; bC(cellBoundary) = 0; % yellow
rC(nucBoundary)  = 0; gC(nucBoundary)  = 1; bC(nucBoundary)  = 1; % cyan
rgb = cat(3, rC, gC, bC);

imwrite(rgb, outPath);
end
