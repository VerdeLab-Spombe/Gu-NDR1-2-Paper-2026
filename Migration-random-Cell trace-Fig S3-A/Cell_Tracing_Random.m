%% Cell migration tracing (manual axes, ASCII-only)
% - Reads X/Y from Excel (header-detect, fallbacks to 3-col/2-col blocks)
% - Removes NaNs, recenters per track, scales px -> um
% - Plots traces with axes at origin (fallback lines), no outer frame
% - Shows numeric tick labels with "nice" spacing
% - Manual axis limits: set xlims/ylims below; use [] for auto.

clear; close all; clc

% ====== CONFIG ======
xlsxFile = 'Scr-shRNA-XY Coordinate-B3-2 -3-3.xlsx';  % change if needed
px2um    = 1.3043;                          % micrometers per pixel
figureID = 1;

% Set manual axis limits in micrometers; use [] for auto
xlims = [-200 200];    % e.g., [-300 500] or []
ylims = [-200 200];    % e.g., [-150 250] or []

% ====== LOAD TRAJECTORIES (pixels, relative to per-track origin) ======
[cellPos_px, origin_px, info] = loadXYTrajectories(xlsxFile); %#ok<NASGU>

% ====== CONVERT TO MICROMETERS ======
cellPos_um = cell(size(cellPos_px));
for ii = 1:numel(cellPos_px)
    P = cellPos_px{ii};
    if ~isempty(P)
        cellPos_um{ii} = P * px2um;
    else
        cellPos_um{ii} = [];
    end
end

% ====== PLOT (units = micrometers) with manual axis limits ======
plotTrajectoriesUM_CustomAxes(cellPos_um, figureID, xlims, ylims);

% ====== METRICS (units = micrometers, on relative traces) ======
metrics = computeMigrationMetrics(cellPos_um);
disp(metrics)
% writetable(metrics, 'cell_migration_metrics_um.csv');  % optional export


% ==================== FUNCTIONS ====================

function [cellPos, origin_px, info] = loadXYTrajectories(xlsxFile)
    % Returns:
    %   cellPos   : cell array of Nx2 doubles (relative traces in pixels)
    %   origin_px : cell array of 1x2 origins (first valid absolute point, pixels)
    % Detects X/Y by headers; falls back to 3-col [X,Y,extra] or 2-col [X,Y] blocks.

    T = readtable(xlsxFile, 'VariableNamingRule','preserve');
    ncol   = size(T,2);
    vnames = T.Properties.VariableNames;   % typically cell array of char

    % Try header-based X/Y pairs
    xyPairs = detectXYPairsByHeader(vnames);

    % Fallback 1: 3-col blocks after first column: [X, Y, extra]
    if isempty(xyPairs)
        if ncol >= 4
            cellNum = floor((ncol - 1)/3);
            tmp = zeros(cellNum,2);
            for k = 1:cellNum
                xidx = k*3;
                yidx = k*3 + 1;
                if yidx <= ncol
                    tmp(k,:) = [xidx, yidx];
                end
            end
            xyPairs = tmp(all(tmp>0,2),:);
        end
    end

    % Fallback 2: 2-col blocks after first column: [X, Y]
    if isempty(xyPairs)
        if ncol >= 3
            cellNum = floor((ncol - 1)/2);
            tmp = zeros(cellNum,2);
            for k = 1:cellNum
                xidx = 1 + (k-1)*2 + 1;  % skip first column (e.g., time)
                yidx = xidx + 1;
                if yidx <= ncol
                    tmp(k,:) = [xidx, yidx];
                end
            end
            xyPairs = tmp(all(tmp>0,2),:);
        end
    end

    cellPos   = {};
    origin_px = {};
    keep = false(size(xyPairs,1),1);

    for i = 1:size(xyPairs,1)
        cols = xyPairs(i,:);
        if any(cols < 1) || any(cols > ncol), continue, end

        Praw = table2array(T(:, cols));   % Nx2 absolute positions (pixels)
        if isempty(Praw), continue, end

        % Remove rows with NaN in either X or Y
        Praw(any(isnan(Praw),2),:) = [];
        if isempty(Praw), continue, end

        % Origin = first valid absolute point (pixels)
        orig = Praw(1,:);
        Ptr  = Praw - orig;               % relative trace (pixels)

        origin_px{end+1} = orig;          %#ok<AGROW>
        cellPos{end+1}   = Ptr;           %#ok<AGROW>
        keep(i) = true;
    end

    info = struct();
    info.VariableNames   = vnames;
    info.DetectedXYPairs = [];
    if any(keep)
        info.DetectedXYPairs = xyPairs(keep,:);
    end
end

function xyPairs = detectXYPairsByHeader(vnames)
    % Char-only version (robust with unique on cellstr)
    % Pairs columns that look like base_X and base_Y (case-insensitive).
    % Also supports X1/Y1 and <base>_x<num>/<base>_y<num>.

    % Ensure vnames is a cell array of char
    if isstring(vnames)
        vnames = cellstr(vnames);
    end

    n = numel(vnames);
    cleaned = cell(n,1);
    for i = 1:n
        s = vnames{i};
        if isstring(s), s = char(s); end
        s = lower(strtrim(regexprep(s, '\s+', ' ')));
        cleaned{i} = s;
    end

    bases = repmat({''}, n, 1);
    isX = false(n,1);
    isY = false(n,1);

    for i = 1:n
        s = cleaned{i};

        % Ignore common non-position columns
        if contains(s, 'time') || contains(s, 'frame') || contains(s, 'id')
            continue
        end

        % case: base + _x or _y (or space/dash)
        tok = regexp(s, '^(.*?)[ _\-]*([xy])$', 'tokens', 'once');
        if ~isempty(tok)
            bases{i} = strtrim(tok{1});
            if strcmp(tok{2}, 'x'), isX(i)=true; else, isY(i)=true; end
            continue
        end

        % case: xNN / yNN
        tok = regexp(s, '^([xy])(\d+)$', 'tokens', 'once');
        if ~isempty(tok)
            bases{i} = tok{2};
            if strcmp(tok{1}, 'x'), isX(i)=true; else, isY(i)=true; end
            continue
        end

        % case: base + _xNN / _yNN
        tok = regexp(s, '^(.*?)[ _\-]*([xy])(\d+)$', 'tokens', 'once');
        if ~isempty(tok)
            bases{i} = sprintf('%s_%s', strtrim(tok{1}), tok{3});
            if strcmp(tok{2}, 'x'), isX(i)=true; else, isY(i)=true; end
            continue
        end
    end

    % Make sure bases is cell array of char, then unique
    bases = cellfun(@char, bases, 'UniformOutput', false);
    [ubases, ~, ib] = unique(bases);

    xyPairs = [];
    for k = 1:numel(ubases)
        if isempty(ubases{k}), continue, end
        idxs = find(ib == k);
        xidx = idxs(isX(idxs));
        yidx = idxs(isY(idxs));
        if numel(xidx) == 1 && numel(yidx) == 1
            xyPairs(end+1,:) = [xidx, yidx]; %#ok<AGROW>
        end
    end

    if ~isempty(xyPairs)
        xyPairs = sortrows(xyPairs, 1);
    end
end

function plotTrajectoriesUM_CustomAxes(cellPos_um, figureID, xlims, ylims)
    % Plot trajectories in micrometers; manual axis limits if provided.
    % xlims and ylims are [min max] in micrometers; use [] for auto.

    figure(figureID); clf
    set(gcf, 'Position', [520, 260, 560, 560]);
    ax = gca; hold(ax, 'on'); axis(ax, 'equal');

    % Remove outer frame, ticks outward
    set(ax, 'Box', 'off', 'TickDir', 'out');

    % Try to place axes at origin (R2020b+). Fallback draws zero lines.
    ok = true;
    try
        set(ax, 'XAxisLocation','origin', 'YAxisLocation','origin');
    catch
        ok = false;
    end
    if ~ok
        try
            xline(0, ':'); yline(0, ':');
        catch
            plot(get(ax,'XLim'), [0 0], ':');
            plot([0 0], get(ax,'YLim'), ':');
        end
    end

    % Gather all points for auto limits
    nonEmpty = ~cellfun(@isempty, cellPos_um);
    allXY = [];
    if any(nonEmpty)
        allXY = cell2mat(cellPos_um(nonEmpty)');
    end

    % Resolve x limits
    if isempty(xlims)
        if ~isempty(allXY)
            xr = range(allXY(:,1));
            pad = 0.05 * max(1, xr);
            xl = [min(allXY(:,1))-pad, max(allXY(:,1))+pad];
        else
            xl = [-100, 100];
        end
    else
        xl = xlims;
    end

    % Resolve y limits
    if isempty(ylims)
        if ~isempty(allXY)
            yr = range(allXY(:,2));
            pad = 0.05 * max(1, yr);
            yl = [min(allXY(:,2))-pad, max(allXY(:,2))+pad];
        else
            yl = [-100, 100];
        end
    else
        yl = ylims;
    end

    % Validate limits
    if ~(numel(xl)==2 && xl(1) < xl(2))
        error('xlims must be [xmin xmax] with xmin < xmax, or [].');
    end
    if ~(numel(yl)==2 && yl(1) < yl(2))
        error('ylims must be [ymin ymax] with ymin < ymax, or [].');
    end

    xlim(xl); ylim(yl);

    % Ticks with numeric labels (nice spacing)
    targetTickCount = 9;  % adjust density
    [xt, xfmt] = niceTicks1D(xl, targetTickCount);
    [yt, yfmt] = niceTicks1D(yl, targetTickCount);
    set(ax, 'XTick', xt, 'YTick', yt);

    % xtickformat/ytickformat may not exist in very old MATLAB; guard them
    try, xtickformat(xfmt); catch, %#ok<CTCH>
        % fallback: leave default tick labels
    end
    try, ytickformat(yfmt); catch, %#ok<CTCH>
    end

    % Plot each track and end point
    for ii = 1:numel(cellPos_um)
        P = cellPos_um{ii};
        if isempty(P), continue, end
        if size(P,1) >= 2
            plot(P(:,1), P(:,2), 'k-');
            plot(P(end,1), P(end,2), 'k.', 'MarkerSize', 18);
        else
            plot(P(1,1), P(1,2), 'k.', 'MarkerSize', 18);
        end
    end

    title('Cell Migration Traces (\mum)');
    xlabel('X (\mum)');
    ylabel('Y (\mum)');

    % Optional: warn if data exceed manual limits
    if ~isempty(allXY)
        if any(allXY(:,1) < xl(1)) || any(allXY(:,1) > xl(2))
            fprintf('Warning: some X data are outside [%g, %g] um.\n', xl(1), xl(2));
        end
        if any(allXY(:,2) < yl(1)) || any(allXY(:,2) > yl(2))
            fprintf('Warning: some Y data are outside [%g, %g] um.\n', yl(1), yl(2));
        end
    end
end

function [ticks, fmt] = niceTicks1D(lims, targetCount)
    % Compute "nice" tick positions and a matching numeric format.
    lo = min(lims); hi = max(lims);
    if ~isfinite(lo) || ~isfinite(hi) || lo == hi
        lo = lo - 1; hi = hi + 1;
    end
    rngv = hi - lo;
    rough = rngv / max(1, targetCount);
    mag = 10.^floor(log10(rough));
    bases = [1, 2, 5, 10];
    step = bases(find(bases.*mag >= rough, 1, 'first'));
    if isempty(step)
        step = 10*mag;
    else
        step = step * mag;
    end

    startTick = ceil(lo/step)*step;
    endTick   = floor(hi/step)*step;
    if startTick > endTick
        startTick = lo; endTick = hi;
    end
    ticks = startTick:step:endTick;

    dec = max(0, -floor(log10(step)));
    dec = min(dec, 3);  % limit decimals
    fmt = sprintf('%%.%df', dec);
end

function M = computeMigrationMetrics(cellPos)
    % Metrics in same units as input (here, micrometers, relative traces).
    % N: time points; NetDisp: net displacement; PathLen: total path length;
    % Meander: NetDisp/PathLen; EndX/EndY: final coordinates.

    n = numel(cellPos);
    N = zeros(n,1);
    NetDisp = zeros(n,1);
    PathLen = zeros(n,1);
    Meander = zeros(n,1);
    EndX = zeros(n,1);
    EndY = zeros(n,1);

    for i = 1:n
        P = cellPos{i};
        if isempty(P)
            N(i) = 0; NetDisp(i) = NaN; PathLen(i) = NaN; Meander(i) = NaN; EndX(i) = NaN; EndY(i) = NaN;
            continue
        end
        N(i) = size(P,1);
        EndX(i) = P(end,1);
        EndY(i) = P(end,2);
        NetDisp(i) = hypot(EndX(i), EndY(i));
        if size(P,1) >= 2
            d = sqrt(sum(diff(P).^2, 2));
            PathLen(i) = sum(d);
            if PathLen(i) > 0
                Meander(i) = NetDisp(i) / PathLen(i);
            else
                Meander(i) = NaN;
            end
        else
            PathLen(i) = 0;
            Meander(i) = NaN;
        end
    end

    M = table((1:n)', N, NetDisp, PathLen, Meander, EndX, EndY, ...
        'VariableNames', {'CellID','N','NetDisp_um','PathLen_um','Meander','EndX_um','EndY_um'});
end
