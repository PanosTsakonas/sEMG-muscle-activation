%% Build LME-ready long table from activations_export.xlsx (robust; skips empty blocks like P5-T5)
infile  = 'activations_export_marginal.xlsx';
outfile = 'activations_long_raw_marginal.csv';

[~, sheets] = xlsfinfo(infile);

% Growable outputs (use numeric + cellstr for max compatibility)
participant_all = [];
trial_all       = [];
sample_all      = [];
muscle_all      = {};   % cellstr
pipeline_all    = {};   % cellstr
activation_all  = [];

skipped = [];  % [participant, trial, reasonCode]

for s = 1:numel(sheets)
    sheetName = sheets{s};
    T = readcell(infile, 'Sheet', sheetName);  % mixed types OK

    r = 1; 
    nrows = size(T,1);

    while r <= nrows
        cell0 = T{r,1};
        if (ischar(cell0) && ~isempty(strtrim(cell0))) || (isstring(cell0) && strlength(cell0)>0)
            cellstr = char(string(cell0));
            if contains(cellstr,'Participant') && contains(cellstr,'Trial')
                % Parse "Participant X - Trial Y"
                tok = regexp(cellstr, 'Participant\s+(\d+)\s*-\s*Trial\s+(\d+)', 'tokens', 'once');
                if isempty(tok), r = r + 1; continue; end
                participant = str2double(tok{1});
                trial       = str2double(tok{2});

                % Header and data block
                hdrRow    = r + 1;
                if hdrRow > nrows
                    skipped = [skipped; participant, trial, 1]; %#ok<AGROW> % no header/data
                    break
                end
                headers   = rowToCellstr(T(hdrRow,:));
                dataStart = hdrRow + 1;
                dataEnd   = dataStart;
                while dataEnd <= nrows
                    rowvals = T(dataEnd,:);
                    if isBlankRow(rowvals), break; end
                    dataEnd = dataEnd + 1;
                end

                % If there are no data rows under the header, skip
                if dataEnd <= dataStart
                    skipped = [skipped; participant, trial, 2]; %#ok<AGROW> % empty block
                    r = dataEnd + 1; 
                    continue
                end

                dataBlock = T(dataStart:dataEnd-1, :);

                % Columns we want
                want = {'Sample','EDC_new','EDC_old','FDS_new','FDS_old','FDP_new','FDP_old'};
                presentMask = ismember(want, headers);
                present = want(presentMask);
                [~, idxMap] = ismember(present, headers);

                % Convert each present column to numeric safely
                S = struct();
                for k = 1:numel(present)
                    S.(present{k}) = toNumericCol(dataBlock(:, idxMap(k)));
                end
                % Ensure Sample exists (1..N)
                if ~isfield(S,'Sample')
                    S.Sample = (1:size(dataBlock,1))';
                end
                Samp = makeIntSafe(S.Sample);

                % Check if any activation columns are present and non-NaN
                actNames = {'EDC_new','EDC_old','FDS_new','FDS_old','FDP_new','FDP_old'};
                haveAny = false; nonNaN = false;
                for c = 1:numel(actNames)
                    if isfield(S, actNames{c})
                        haveAny = true;
                        if any(~isnan(S.(actNames{c})))
                            nonNaN = true;
                            break
                        end
                    end
                end
                if ~haveAny
                    skipped = [skipped; participant, trial, 3]; %#ok<AGROW> % no activation columns
                    r = dataEnd + 1; 
                    continue
                end
                if ~nonNaN
                    skipped = [skipped; participant, trial, 4]; %#ok<AGROW> % all NaN activations
                    r = dataEnd + 1; 
                    continue
                end

                % Append long rows for each muscle/pipeline
                musList = {'EDC','FDS','FDP'};
                for mu = 1:numel(musList)
                    for flagNew = [true false]
                        % Build column name suffix and pipeline label
                        if flagNew
                            suffix  = 'new';
                            plabel  = 'New';
                        else
                            suffix  = 'old';
                            plabel  = 'Old';
                        end
                        colName = [musList{mu} '_' suffix];

                        if isfield(S, colName)
                            vals  = S.(colName);
                            valid = ~isnan(vals);
                            n     = nnz(valid);
                            if n == 0, continue; end

                            participant_all = [participant_all; repmat(participant, n, 1)]; %#ok<AGROW>
                            trial_all       = [trial_all;       repmat(trial,       n, 1)]; %#ok<AGROW>
                            sample_all      = [sample_all;      Samp(valid)]; %#ok<AGROW>
                            muscle_all      = [muscle_all;      repmat({musList{mu}}, n, 1)]; %#ok<AGROW>
                            pipeline_all    = [pipeline_all;    repmat({plabel}, n, 1)]; %#ok<AGROW>
                            activation_all  = [activation_all;  vals(valid)]; %#ok<AGROW>
                        end
                    end
                end

                r = dataEnd + 1; 
                continue
            end
        end
        r = r + 1;
    end
end

% Assemble table
longTbl = table( ...
    participant_all, trial_all, sample_all, muscle_all, pipeline_all, activation_all, ...
    'VariableNames', {'participant','trial','sample','muscle','pipeline','activation'});

% Normalize time within each (participant, trial) to [0,1]
[~, ~, grpIdx] = unique(longTbl(:,{'participant','trial'}), 'rows');
minS = accumarray(grpIdx, double(longTbl.sample), [], @min);
maxS = accumarray(grpIdx, double(longTbl.sample), [], @max);
time_std = zeros(height(longTbl),1);
for g = 1:numel(minS)
    msk = (grpIdx == g);
    den = max(1, maxS(g) - minS(g));
    time_std(msk) = (double(longTbl.sample(msk)) - minS(g)) ./ den;
end
longTbl.time_std = time_std;
longTbl = movevars(longTbl, 'time_std', 'After', 'sample');

% Write CSV
writetable(longTbl, outfile);
fprintf('Wrote %s with %d rows.\n', outfile, height(longTbl));

% Report skipped blocks (e.g., P5–T5)
if ~isempty(skipped)
    reason = {'no header/data','empty block','no activation columns','all NaN activations'};
    for i = 1:size(skipped,1)
        fprintf('Skipped P%d Trial %d (%s)\n', skipped(i,1), skipped(i,2), reason{skipped(i,3)});
    end
end

%% ---------- Local helper functions ----------
function vs = rowToCellstr(rowCells)
    % Convert a mixed-type header row to a cell array of char
    n = numel(rowCells);
    vs = cell(1,n);
    for i = 1:n
        xi = rowCells{i};
        if isstring(xi)
            vs{i} = char(xi);
        elseif ischar(xi)
            vs{i} = xi;
        elseif isnumeric(xi) && ~isempty(xi) && ~isnan(xi)
            vs{i} = char(string(xi));
        else
            vs{i} = '';
        end
    end
end

function tf = isBlankRow(rowCells)
    % True if row is empty/blank/NaN across columns
    tf = true;
    for i = 1:numel(rowCells)
        xi = rowCells{i};
        if isstring(xi)
            if strlength(strtrim(xi)) > 0, tf = false; return; end
        elseif ischar(xi)
            if ~isempty(strtrim(xi)), tf = false; return; end
        elseif isnumeric(xi)
            if ~isempty(xi) && ~all(isnan(xi)), tf = false; return; end
        end
    end
end

function v = toNumericCol(c)
    % Safe numeric conversion for a mixed-type cell column
    n = size(c,1);
    v = nan(n,1);
    for i = 1:n
        x = c{i};
        if isnumeric(x)
            if ~isempty(x), v(i) = x; end
        elseif isstring(x) || ischar(x)
            t = str2double(char(string(x)));
            if ~isnan(t), v(i) = t; end
        else
            % leave NaN
        end
    end
end

function v = makeIntSafe(v)
    % Convert to integer vector safely (NaN already handled upstream)
    v = round(double(v));
    v(~isfinite(v)) = 0;
end
