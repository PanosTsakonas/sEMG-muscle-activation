
%% LME (simpler): log(activation) ~ pipeline + muscle + (1|participant)
% Per-muscle models:   log(activation) ~ pipeline + (1|participant)
%
% Outputs are saved with "_simple" suffix for side-by-side comparison with your
% previous full-model exports.
%
% Files written:
%   - lme_results_summary_simple.csv
%   - lme_fixed_effects_overall_simple.csv
%   - lme_cov_params_overall_simple.csv (or a *_WARNING.txt fallback)
%   - lme_random_effects_overall_simple.csv (or .txt fallback)
%   - lme_fit_stats_overall_simple.csv
%   - lme_summary_overall_simple.txt
%   - ... and the same 5 files for EDC, FDS, FDP
close all
clear all
infile  = 'activations_long_raw_marginal.csv';
summary_outcsv  = 'LME_Results/lme_results_summary_simple.csv';

% ========= Settings =========
USE_BINNING      = true;   % keep true for comparability with your previous run
NBINS_PER_TRIAL  = 20;     % bins per (participant×trial×muscle×pipeline)
EPS              = 1e-6;   % clamp for log-transform

% ========= Load & prep =========
T = readtable(infile);
req = {'participant','trial','sample','time_std','muscle','pipeline','activation'};
assert(all(ismember(req, T.Properties.VariableNames)), 'CSV missing required columns.');

% Clamp zeros/negatives, log-transform
T.activation = max(T.activation, EPS);
Temp     = log(T.activation);
% Example: Winsorisation of data vector 'dataVec' at 5% tails
p = 0.05;

% Calculate lower and upper percentile bounds
lower_bound = quantile(Temp, p);
upper_bound = quantile(Temp, 1 - p);

% Winsorise: replace values beyond bounds with bound values
winsorisedData = Temp;
winsorisedData(winsorisedData < lower_bound) = lower_bound;
winsorisedData(winsorisedData > upper_bound) = upper_bound;
T.logAct=(winsorisedData);

% Categoricals
T.participant = categorical(T.participant);
T.trial       = categorical(T.trial);
T.muscle      = categorical(T.muscle, {'EDC','FDS','FDP'});
T.pipeline    = categorical(T.pipeline, {'Old','New'});  % Old = baseline

% ========= Optional binning (same as before for fair comparison) =========
if USE_BINNING
    fprintf('Binning each (participant×trial×muscle×pipeline) to %d bins ...\n', NBINS_PER_TRIAL);
    T = bin_by_time_fast(T, NBINS_PER_TRIAL);
end

% Sanity checks
fprintf('N rows: %d | Participants: %d\n', height(T), numel(categories(T.participant)));
disp('Pipelines:'), disp(categories(T.pipeline)')
disp('Muscles  :'), disp(categories(T.muscle)')

% % ========= Overall simpler model =========
% % logAct ~ pipeline + muscle + (1|participant)
% formula_overall = 'logAct ~ pipeline + muscle + trial + (1|participant)';
% lme_overall = fitlme(T, formula_overall, 'FitMethod','REML')
% 
% 
% [ratio, rlow, rhigh, b_log, ci_log] = pipeline_ratio_from_model(lme_overall, 'pipeline_New');
% 
% fprintf('\n=== OVERALL (simple) New vs Old ===\n');
% fprintf('GMR (New/Old) = %.3f (95%% CI %.3f–%.3f)\n', ratio, rlow, rhigh);
% 
 rows = {};
% rows = [rows; {"Overall", b_log, ci_log(1), ci_log(2), ratio, rlow, rhigh}];
% 
% % Export full details for overall (with "_simple" suffix)
% export_full_model(lme_overall, 'overall_simple');

% ========= Per-muscle simpler models =========
muscles = categories(T.muscle);
for m = 1:numel(muscles)
    mus = muscles{m};
    subT = T(T.muscle == mus,:);
    % logAct ~ pipeline + (1|participant)
    formula_mus = 'logAct ~ pipeline + trial + (1|participant)';
    disp(muscles{m});
    lme_m = fitlme(subT, formula_mus, 'FitMethod','REML')
    %lme_m = fitglme(subT, formula_mus, 'Distribution', 'Normal','Link','log');
    % Check residual normality via histogram and qqplot
    res = residuals(lme_m);
    % Fit normal distribution
    pd = fitdist(res, 'Normal');
    figure;
    histogram(res,'Normalization','pdf'); hold on;
%     [h,p]=swtest(res);
%     legend("KS test p-value: "+p);
    x_values = linspace(min(res), max(res), 100);
    plot(x_values, pdf(pd, x_values), 'r-', 'LineWidth', 2);
    legend("Residuals","Normal distribution");
    title("Residual Histogram of "+mus+" muscle");
    saveas(gcf,"LME_Results/Marginal/histogram_residuals_muscle_"+mus+".png");
    figure()
     qqplot(res);
     title("Normal Q-Q Plot of "+mus+" muscle");
    saveas(gcf,"LME_Results/Marginal/qq_plot_residuals_muscle_"+mus+".png");

    % Use Holm-Bonferroni correction for multiple testing on fixed effects p-values
    fixedEffectsTbl = lme_m.Coefficients;
    pvals = fixedEffectsTbl.pValue;
    [p_sorted,idx] = sort(pvals);
    N = length(pvals);
    adj_pvals = zeros(size(pvals));
    for i = 1:N
        adj_pvals(idx(i)) = min((N - i + 1)*p_sorted(i),1);
    end
    fixedEffectsTbl.adj_pValue = adj_pvals;
    disp(fixedEffectsTbl);
    if ~istable(fixedEffectsTbl)
        fixedEffectsTbl = dataset2table(fixedEffectsTbl);
    end
    % Save adjusted table for reviewer
    writetable(fixedEffectsTbl, 'LME_Results/Marginal/lme_fixed_effects_adj_pvalues.csv');

    [ratio, rlow, rhigh, b_log, ci_log] = pipeline_ratio_from_model(lme_m, 'pipeline_New');

    fprintf('\n=== %s (simple) New vs Old ===\n', mus);
    fprintf('GMR (New/Old) = %.3f (95%% CI %.3f–%.3f)\n', ratio, rlow, rhigh);

    rows = [rows; {mus, b_log, ci_log(1), ci_log(2), ratio, rlow, rhigh}];

    % Export full details per muscle (with "_simple" suffix)
    export_full_model(lme_m, [lower(mus) '_simple']);
    clear lme_m
end

% ========= Save summary (ratios) =========
VarNames = {'Model','Beta_log','CI_low_log','CI_high_log','GMR','CI_low_ratio','CI_high_ratio'};
Summary = cell2table(rows,'VariableNames',VarNames);
writetable(Summary, summary_outcsv);
fprintf('\nSaved simple-model ratios to %s\n', summary_outcsv);



%% ---------- Helper: extract New vs Old ratio robustly ----------
function [ratio, rlow, rhigh, b_log, ci_log] = pipeline_ratio_from_model(lme, termName)
    coefTbl = ensureTable(lme.Coefficients);    % robust conversion to table
    ciMat   = coefCI(lme);                      % same row order as Coefficients

    idx = [];
    if ismember('Name', coefTbl.Properties.VariableNames)
        idx = find(strcmp(coefTbl.Name, termName), 1);
    end
    if isempty(idx) && ~isempty(coefTbl.Properties.RowNames)
        idx = find(strcmp(coefTbl.Properties.RowNames, termName), 1);
    end
    if isempty(idx)
        for k = 1:width(coefTbl)                 % last resort: any string column
            col = coefTbl{:,k};
            if iscellstr(col) || isstring(col)
                f = find(strcmp(cellstr(col), termName), 1);
                if ~isempty(f), idx = f; break; end
            end
        end
    end
    assert(~isempty(idx), 'Could not find term "%s" in model coefficients.', termName);

    est = tryGetNumeric(coefTbl, {'Estimate','Est','Beta','Value'});
    b_log = est(idx);
    ci_row = ciMat(idx, :);
    ci_log = [ci_row(1), ci_row(2)];

    ratio  = exp(b_log);
    rlow   = exp(ci_log(1));
    rhigh  = exp(ci_log(2));
end

%% ---------- Helper: export full model details (robust across releases) ----------
function export_full_model(lme, tag)
    % 1) Fixed effects
    fixedCSV = sprintf('LME_Results/lme_fixed_effects_%s.csv', tag);
    try
        coefTbl = ensureTable(lme.Coefficients);
        writetable(coefTbl, fixedCSV);
    catch ME
        warnTxt = sprintf('WARNING: Failed to export fixed effects for "%s": %s', tag, ME.message);
        fprintf('%s\n', warnTxt);
        fid = fopen(strrep(fixedCSV,'.csv','_WARNING.txt'),'w'); fprintf(fid,'%s\n',warnTxt); fclose(fid);
    end

    % 2) Covariance parameters (random-effects variances & residual)
    covCSV = sprintf('lme_cov_params_%s.csv', tag);
    covTbl = safeCovParams(lme);  % handles property OR method OR returns []
    if ~isempty(covTbl)
        try
            writetable(covTbl, covCSV);
        catch
            fid = fopen(strrep(covCSV,'.csv','.txt'),'w');
            fprintf(fid,'%s\n', evalc('disp(covTbl)')); 
            fclose(fid);
        end
    else
        warnTxt = sprintf('WARNING: Covariance parameters not available on this MATLAB release for "%s".', tag);
        fprintf('%s\n', warnTxt);
        fid = fopen(strrep(covCSV,'.csv','_WARNING.txt'),'w'); fprintf(fid,'%s\n',warnTxt); fclose(fid);
    end

    % 3) Random effects (participant BLUPs)
    reCSV = sprintf('lme_random_effects_%s.csv', tag);
    try
        RE = randomEffects(lme);              % may be object/dataset/table
        try
            REtbl = ensureTable(RE);
            writetable(REtbl, reCSV);
        catch
            reTXT = strrep(reCSV,'.csv','.txt');     % fallback: text dump
            fid = fopen(reTXT,'w'); fprintf(fid,'%s\n', evalc('disp(RE)')); fclose(fid);
            fprintf('Saved random effects as text (fallback) for "%s": %s\n', tag, reTXT);
        end
    catch
        warnTxt = sprintf('WARNING: randomEffects() not available on this MATLAB release for "%s".', tag);
        fprintf('%s\n', warnTxt);
        fid = fopen(strrep(reCSV,'.csv','_WARNING.txt'),'w'); fprintf(fid,'%s\n',warnTxt); fclose(fid);
    end

    % 4) Fit statistics
    fitCSV = sprintf('lme_fit_stats_%s.csv', tag);
    fitStats = table(lme.LogLikelihood, lme.ModelCriterion.AIC, lme.ModelCriterion.BIC, ...
                     lme.ModelCriterion.Deviance, ...
                     'VariableNames', {'LogLikelihood','AIC','BIC','Deviance'});
    writetable(fitStats, fitCSV);

    % 5) Text summary (ANOVA + printed model)
    txtPath = sprintf('lme_summary_%s.txt', tag);
    s = evalc('disp(lme)');        % printed model summary
    try
        a = evalc('anova(lme)');   % Wald tests for fixed effects
        s = sprintf('%s\n\nANOVA (Wald tests):\n%s', s, a);
    catch
        % anova may be missing; ignore
    end
    fid = fopen(txtPath, 'w'); fwrite(fid, s); fclose(fid);

    fprintf('Exported details for "%s":\n  %s\n  %s\n  %s\n  %s\n  %s\n', ...
        tag, fixedCSV, covCSV, reCSV, fitCSV, txtPath);
end

%% ---------- Helper: robust covariance-parameters to table ----------
function tbl = safeCovParams(lme)
    % Try property (newer releases)
    tbl = [];
    try
        CP = lme.CovarianceParameters; %#ok<NASGU>
        tbl = ensureTable(lme.CovarianceParameters);
        return
    catch
    end
    % Try method form (older releases often provide covarianceParameters(lme))
    try
        try
            [Psi, Names] = covarianceParameters(lme); %#ok<ASGLU>
        catch
            Psi = covarianceParameters(lme); Names = [];
        end
        if iscell(Psi)
            rows = {};
            for g = 1:numel(Psi)
                Pg = Psi{g};
                if isempty(Pg), continue; end
                vals = Pg(:);
                for k = 1:numel(vals)
                    rows = [rows; {g, k, vals(k)}]; %#ok<AGROW>
                end
            end
            tbl = cell2table(rows, 'VariableNames', {'GroupIndex','ParamIndex','Estimate'});
        elseif isnumeric(Psi)
            tbl = table((1:numel(Psi))', Psi(:), 'VariableNames', {'ParamIndex','Estimate'});
        else
            tbl = ensureTable(Psi);
        end
        return
    catch
        tbl = [];
        return
    end
end

%% ---------- Helper: robustly convert tables/datasets/objects to table ----------
function tbl = ensureTable(x)
    if istable(x), tbl = x; return, end
    if exist('dataset2table','file') == 2
        try, tbl = dataset2table(x); return; catch, end
    end
    try, tbl = struct2table(x); return; catch, end
    try
        fields = {'Name','Estimate','SE','tStat','DF','pValue','Group','Level'};
        present = intersect(fields, fieldnames(x));
        if ~isempty(present)
            S = struct();
            for i = 1:numel(present), S.(present{i}) = x.(present{i}); end
            tbl = struct2table(S); return
        end
    catch
    end
    error('Could not convert object to table for export on this MATLAB release.');
end

%% ---------- Helper: safely pull a numeric column by candidate names ----------
function col = tryGetNumeric(tbl, names)
    for k = 1:numel(names)
        if ismember(names{k}, tbl.Properties.VariableNames)
            v = tbl.(names{k});
            if isnumeric(v), col = v; return, end
        end
    end
    for k = 1:width(tbl)
        v = tbl{:,k};
        if isnumeric(v), col = v; return, end
    end
    error('No numeric estimate column found in coefficient table.');
end

%% ---------- Helper: fast binning via discretize + groupsummary ----------
function Tout = bin_by_time_fast(T, nbins)
    % Bin index per (participant, trial, muscle, pipeline)
    edges = linspace(0, 1, nbins+1);
    G = findgroups(T.participant, T.trial, T.muscle, T.pipeline);

    bin = NaN(height(T),1);
    ug = unique(G);
    for i = 1:numel(ug)
        mask = (G == ug(i));
        ts = T.time_std(mask);
        ts(ts == 1) = 1 - eps;               % keep right-edge in last bin
        bin(mask) = discretize(ts, edges);
    end

    T.bin = bin;
    grpvars = {'participant','trial','muscle','pipeline','bin'};
    Agg = groupsummary(T, grpvars, "mean", {'time_std','activation','logAct'});

    Agg.Properties.VariableNames{'mean_time_std'} = 'time_std';
    Agg.Properties.VariableNames{'mean_activation'} = 'activation';
    Agg.Properties.VariableNames{'mean_logAct'}    = 'logAct';

    Tout = Agg(:, {'participant','trial','muscle','pipeline','time_std','activation','logAct'});
end
