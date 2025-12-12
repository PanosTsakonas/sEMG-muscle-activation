clc;
close all;


fsEMG=1000;
fn=fsEMG/2;

file=*************Set the path to the sEMG Data\*********************";

%Set 1 for zero initial condition in the ODE solver or to [] for non zero
Zero=1;

if isempty(Zero)==0

    ff=file+"Plots\Activation_zero_init_VMD\";
   
    ffr=file+"Plots\Ratio_zero_init_VMD\";
    
else
    ff=file+"Plots\Activation_nonzero_init\";
    ffr=file+"Plots\Ratio_nonzero_init\";
end

%This checks whether the folders exsit that will have all the plots
File_Create(ff);
File_Create(ffr);

ff_Sig=file+"Plots\Signal_VMD\";
File_Create(ff_Sig);

ff_ODE=file+"Plots\Activation ODE zero init both VMD\";
ff_ODE_R=file+"Plots\Ratio of Activation ODE zero init both VMD\";
f_box="LME_Results\intra_subject_Wilcoxon_New\";
f_inter=file+"Plots\BoxPlots\inter_subject VMD\";
File_Create(f_box);
File_Create(ff_ODE);
File_Create(ff_ODE_R);
File_Create(f_inter);

%Set to a number to allow the HHT plot of sEMGs
HHT_View=1;

if isempty(HHT_View)==0
    ff_HHT="C:\Users\panos\Desktop\Plots\HHT_Marginal\";
    File_Create(ff_HHT);
end

%Set the filter



[bf,af]=butter(4,30/fn,'low');


mov_avg_window = 0.05; % Moving average window size (seconds) taken from 
% "An EMG-to-Force Processing Approach to Estimating Knee Muscle Forces during Adult, 
% Self-Selected Speed Gait, where they apply a moving average over 50
% samples

%Set MVC

MVC=zeros(10,3);


for i=1:10

        D=importfile(file+"sEMG_P"+i+".xlsx",6);
        MVC(i,:)=D.data;
end

c=1;
for i=1:10

    for j=1:5
         t=[];
        if (i==5 && j==5)
       edc{c}=[];
       fds{c}=[];
       fdp{c}=[];
       sedc{c}=[];
       sfds{c}=[];
       sfdp{c}=[];
       R_EDC{c}=[];
       R_FDP{c}=[];
       R_FDS{c}=[];
        c=c+1;
        else

        D=importfile(file+"sEMG_P"+i+".xlsx",j);

        EDC=double(D.data(:,1));
        FDS=double(D.data(:,2));
        FDP=double(D.data(:,3));

        %HHT
        [hsEDC,f_EDC,tsEDC]=hht(emd(EDC),fsEMG);
        [hsFDS,f_FDS,tsFDS]=hht(emd(FDS),fsEMG);
        [hsFDP,f_FDP,tsFDP]=hht(emd(FDP),fsEMG);

        [bEDC,aEDC,fEDC]=HHT_Filter_Cut_off(hsEDC,tsEDC,fn,f_EDC,HHT_View,ff_HHT,i,j,"EDC",EDC);
        [bFDP,aFDP,fFDP]=HHT_Filter_Cut_off(hsFDS,tsFDS,fn,f_FDP,HHT_View,ff_HHT,i,j,"FDS",FDS);
        [bFDS,aFDS,fFDS]=HHT_Filter_Cut_off(hsFDP,tsFDP,fn,f_FDS,HHT_View,ff_HHT,i,j,"FDP",FDP);
    
    
        for m=1:length(FDP)

            if m==1
                t(m)=0;
            else
                t(m)=t(m-1)+1/fsEMG;
            end
        end

        
        EDCf=abs(filtfilt(bEDC,aEDC,EDC))./MVC(i,1);
        FDPf=abs(filtfilt(bFDP,aFDP,FDP))./MVC(i,3);
        FDSf=abs(filtfilt(bFDS,aFDS,FDS))./MVC(i,2);
        
        %Get the envelopes
        EDC_env=upper_envelope(EDCf,t);
        FDP_env=upper_envelope(FDPf,t);
        FDS_env=upper_envelope(FDSf,t);

        %Solve the ODE from Blana et al.

        if (isempty(Zero)==0)

            [~,a_edc_ode_env]=ode45(@(t,a) activation(t,a,EDC_env),t,0);
            [~,a_fdp_ode_env]=ode45(@(t,a) activation(t,a,FDP_env),t,0);
            [~,a_fds_ode_env]=ode45(@(t,a) activation(t,a,FDS_env),t,0);
        else
            [~,a_edc_ode_env]=ode45(@(t,a) activation(t,a,EDC_env),t,EDCf(1));
            [~,a_fdp_ode_env]=ode45(@(t,a) activation(t,a,FDP_env),t,FDPf(1));
            [~,a_fds_ode_env]=ode45(@(t,a) activation(t,a,FDS_env),t,FDSf(1));
        end


        %Based on the A fast implementation for EMG signal linear envelope computation

        mov_avg_samples = round(mov_avg_window * fsEMG);
        smooth_EDC = filtfilt(bf,af,movmean(EDCf, mov_avg_samples));
        smooth_FDS = filtfilt(bf,af,movmean(FDSf, mov_avg_samples));
        smooth_FDP = filtfilt(bf,af,movmean(FDPf, mov_avg_samples));

        %Set Spline and solve ODE

        pedc=spline(t,smooth_EDC);
        pfdp=spline(t,smooth_FDP);
        pfds=spline(t,smooth_FDS);

        ppedc=mkpp(pedc.breaks,pedc.coefs);
        ppfdp=mkpp(pfdp.breaks,pfdp.coefs);
        ppfds=mkpp(pfds.breaks,pfds.coefs);

        %Solve the ODE from Blana et al.

        [~,a_edc_ode]=ode45(@(t,a) activation(t,a,ppedc),t,0);
        [~,a_fdp_ode]=ode45(@(t,a) activation(t,a,ppfdp),t,0);
        [~,a_fds_ode]=ode45(@(t,a) activation(t,a,ppfds),t,0);

        %Display outputs

        % Print("EDC",a_edc_ode,a_edc_ode_env,t,ff_ODE,i,j);
        % Print("FDP",a_fdp_ode,a_fdp_ode_env,t,ff_ODE,i,j);
        % Print("FDS",a_fds_ode,a_fds_ode_env,t,ff_ODE,i,j);
        % 
        % PrintR("EDC",a_edc_ode,a_edc_ode_env,t,ff_ODE_R,i,j);
        % PrintR("FDP",a_fdp_ode,a_fdp_ode_env,t,ff_ODE_R,i,j);
        % PrintR("FDS",a_fds_ode,a_fds_ode_env,t,ff_ODE_R,i,j);

            sedc{c}=smooth_EDC.*100;
            sfds{c}=smooth_FDS.*100;
            sfdp{c}=smooth_FDP.*100;
            edc{c}=a_edc_ode_env.*100;
            fds{c}=a_fds_ode_env.*100;
            fdp{c}=a_fdp_ode_env.*100;
        
            R_EDC{c}=a_edc_ode_env./smooth_EDC;
            R_FDS{c}=a_fds_ode_env./smooth_FDS;
            R_FDP{c}=a_fdp_ode_env./smooth_FDP;
            c=c+1;
        

        % figure()
        % plot(t,smooth_EDC,t,ppval(EDC_env,t),'--','LineWidth',1)
        % hold on;
        % plot(t,EDCf);
        % hold off;
        % title("EDC signal");
        % legend("Existing pipeline","New pipeline","EMG signal");
        % xlabel("Time (s)");
        % ylabel("Normalised EMG signal");
        % saveas(gcf,ff_Sig+"EDC_Signal_Par_"+i+"_Cyl_"+j+".png");
        % 
        % figure()
        % plot(t,smooth_FDS,t,ppval(FDS_env,t),'--','LineWidth',1)
        % hold on;
        % plot(t,FDSf);
        % hold off;
        % title("FDS signal");
        % legend("Existing pipeline","New pipeline","EMG signal");
        % xlabel("Time (s)");
        % ylabel("Normalised EMG signal");
        % saveas(gcf,ff_Sig+"FDS_Signal_Par_"+i+"_Cyl_"+j+".png");
        % 
        % figure()
        % plot(t,smooth_FDP,t,ppval(FDP_env,t),'--','LineWidth',1)
        % hold on;
        % plot(t,FDPf);
        % hold off;
        % title("FDP signal");
        % legend("Existing pipeline","New pipeline","EMG signal");
        % xlabel("Time (s)");
        % ylabel("Normalised EMG signal");
        % saveas(gcf,ff_Sig+"FDP_Signal_Par_"+i+"_Cyl_"+j+".png");
        % 
        % figure()
        % plot(t,smooth_EDC.*100,t,a_edc_ode_env.*100)
        % title("EDC muscle activation");
        % legend("Windowing pipeline","New pipeline");
        % xlabel("Time (s)");
        % ylabel("Muscle activation (%)");
        % saveas(gcf,ff+"EDC_Activation_Par_"+i+"_Cyl_"+j+".png");
        % 
        % figure()
        % plot(t,smooth_FDP.*100,t,a_fdp_ode_env.*100)
        % title("FDP muscle activation");
        % legend("Windowing pipeline","New pipeline");
        % xlabel("Time (s)");
        % ylabel("Muscle activation (%)");
        % saveas(gcf,ff+"FDP_Activation_Par_"+i+"_Cyl_"+j+".png");
        % 
        % figure()
        % plot(t,smooth_FDS.*100,t,a_fds_ode_env.*100)
        % title("FDS muscle activation");
        % legend("Windowing pipeline","New pipeline");
        % xlabel("Time (s)");
        % ylabel("Muscle activation (%)");
        % saveas(gcf,ff+"FDS_Activation_Par_"+i+"_Cyl_"+j+".png");
        % 
        %  figure()
        % plot(t,a_edc_ode_env./smooth_EDC)
        % title("Ratio between new and existing pipelines for EDC muscle");
        % xlabel("Time (s)");
        % ylabel("Ratio of muscle activations");
        % saveas(gcf,ffr+"EDC_Ratio_Par_"+i+"_Cyl_"+j+".png");
        % 
        % figure()
        % plot(t,a_fdp_ode_env./smooth_FDP)
        % title("Ratio between new and existing pipelines for FDP muscle");
        % xlabel("Time (s)");
        % ylabel("Ratio of muscle activations");
        % saveas(gcf,ffr+"FDP_Ratio_Par_"+i+"_Cyl_"+j+".png");
        % 
        % figure()
        % plot(t,a_fds_ode_env./smooth_FDS)
        % title("Ratio between new and existing pipelines for FDS muscle");
        % xlabel("Time (s)");
        % ylabel("Ratio of muscle activations");
        % saveas(gcf,ffr+"FDS_Ratio_Par_"+i+"_Cyl_"+j+".png");

        end
        clear t EDC FDP FDS

    end

end

%%
%Intra Subject boxplots
EDC_MED = cell(10,1); 
FDS_MED = cell(10,1); 
FDP_MED = cell(10,1); 

for i=1:10
[M_EDC,sd_EDC,SM_EDC,sd_SEDC,MED_EDC]=intra_Boxplot(edc,sedc,i,"EDC",f_box);
[M_FDS,sd_FDS,SM_FDS,sd_SFDS,MED_FDS]=intra_Boxplot(fds,sfds,i,"FDS",f_box);
[M_FDP,sd_FDP,SM_FDP,sd_SFDP,MED_FDP]=intra_Boxplot(fdp,sfdp,i,"FDP",f_box);


% store participant-level of all_diffs
    EDC_MED{i} = MED_EDC;  
    FDS_MED{i} = MED_FDS;  
    FDP_MED{i} = MED_FDP;  
end

%% Format each cell as "median (IQR)", show "NA" if missing
fmtCell = @(m,i) ...
    ( ...
        (isnan(m) || isnan(i)) * 0 + 0 ... % dummy to satisfy syntax
    ); % placeholder (ignored below)

% Replace with a normal function handle that uses if/else logic:
fmtCell = @(m,i) ...
    ( ...
        (isnan(m) || isnan(i)) ...
        && "NA" ...
        || string(sprintf('%.3f (%.3f)', m, i)) ...
    );

% 10x1
Participant = (1:10)';

% helper to get per-trial scalar differences and pad to 5 with NaN
perTrial = @(M) [ median(M, 1, 'omitnan'), nan(1, max(0, 5 - size(M,2))) ];

EDC_T = nan(10,5);
FDS_T = nan(10,5);
FDP_T = nan(10,5);

for ii = 1:10
    % Each cell may be n x 4 (subject 5) or n x 5 (others). Collapse columns (trials).
    if ~isempty(EDC_MED{ii}), EDC_T(ii,:) = perTrial(EDC_MED{ii}); end
    if ~isempty(FDS_MED{ii}), FDS_T(ii,:) = perTrial(FDS_MED{ii}); end
    if ~isempty(FDP_MED{ii}), FDP_T(ii,:) = perTrial(FDP_MED{ii}); end
end

% Round to 3 decimal places
EDC_T = round(EDC_T, 3);
FDS_T = round(FDS_T, 3);
FDP_T = round(FDP_T, 3);

% Build table with columns: EDC_T1..T5, FDS_T1..T5, FDP_T1..T5
varNames = [{'Participant'}, ...
            arrayfun(@(k) sprintf('EDC_T%d',k), 1:5, 'UniformOutput', false), ...
            arrayfun(@(k) sprintf('FDS_T%d',k), 1:5, 'UniformOutput', false), ...
            arrayfun(@(k) sprintf('FDP_T%d',k), 1:5, 'UniformOutput', false)];

T = table(Participant, ...
    EDC_T(:,1), EDC_T(:,2), EDC_T(:,3), EDC_T(:,4), EDC_T(:,5), ...
    FDS_T(:,1), FDS_T(:,2), FDS_T(:,3), FDS_T(:,4), FDS_T(:,5), ...
    FDP_T(:,1), FDP_T(:,2), FDP_T(:,3), FDP_T(:,4), FDP_T(:,5), ...
    'VariableNames', varNames);

% === Export for the appendix ===
out_csv  = fullfile(f_box, 'AppendixA_median_IQR_differences_Marginal.csv');
writetable(T, out_csv);


%% Create activations_export excel

outfile = "activations_export_marginal.xlsx";

nParticipants = 10;
nTrialsPerP   = 5;

for p = 1:nParticipants
    sheetName = sprintf('P%d', p);
    C = {};          % cell buffer for this sheet
    row = 1;

    for t = 1:nTrialsPerP
        % Linear index into your cell arrays
        c = (p-1)*nTrialsPerP + t;
        if c > numel(edc)
            continue;   % safety
        end

        % Extract data for this participant/trial
        edc_t  = edc{c};
        fds_t  = fds{c};
        fdp_t  = fdp{c};
        sedc_t = sedc{c};
        sfds_t = sfds{c};
        sfdp_t = sfdp{c};

        % Make sure they're column vectors
        if ~isempty(edc_t),  edc_t  = edc_t(:);  end
        if ~isempty(fds_t),  fds_t  = fds_t(:);  end
        if ~isempty(fdp_t),  fdp_t  = fdp_t(:);  end
        if ~isempty(sedc_t), sedc_t = sedc_t(:); end
        if ~isempty(sfds_t), sfds_t = sfds_t(:); end
        if ~isempty(sfdp_t), sfdp_t = sfdp_t(:); end

        % If this trial has absolutely no data, skip (e.g. P5–T5)
        if isempty(edc_t)  && isempty(fds_t)  && isempty(fdp_t) && ...
           isempty(sedc_t) && isempty(sfds_t) && isempty(sfdp_t)
            fprintf('Skipping Participant %d Trial %d (no data)\n', p, t);
            continue;
        end

        % Determine length L of this trial (assume all series same length;
        % to be safe, use min of non-empty lengths)
        lengths = [numel(edc_t), numel(sedc_t), ...
                   numel(fds_t), numel(sfds_t), ...
                   numel(fdp_t), numel(sfdp_t)];
        L = min(lengths(lengths > 0));
        if isempty(L) || L == 0
            fprintf('Skipping Participant %d Trial %d (degenerate lengths)\n', p, t);
            continue;
        end

        % ---------- Block header: "Participant X - Trial Y" ----------
        C{row,1} = sprintf('Participant %d - Trial %d', p, t);
        row = row + 1;

        % ---------- Column headers ----------
        C(row,1) = {'Sample'};
        C(row,2) = {'EDC_new'};
        C(row,3) = {'EDC_old'};
        C(row,4) = {'FDS_new'};
        C(row,5) = {'FDS_old'};
        C(row,6) = {'FDP_new'};
        C(row,7) = {'FDP_old'};
        row = row + 1;

        % ---------- Data rows (no NaN padding, like your original) ----------
        for i = 1:L
            C{row,1} = i;            % Sample index
            C{row,2} = edc_t(i);     % new pipeline
            C{row,3} = sedc_t(i);    % old pipeline
            C{row,4} = fds_t(i);
            C{row,5} = sfds_t(i);
            C{row,6} = fdp_t(i);
            C{row,7} = sfdp_t(i);
            row = row + 1;
        end

        % ---------- Blank row between trial blocks ----------
        C(row,1) = {[]};
        row = row + 1;
    end

    % Only write sheet if it has any content
    if ~isempty(C)
        writecell(C, outfile, 'Sheet', sheetName);
        fprintf('Wrote sheet %s with %d rows.\n', sheetName, row-1);
    end
end

function env=upper_envelope(Signal,t)

    %Create the upper envelope of the signal
        
    %Get the magnitude of the Hilbert transform

    Sig=abs(hilbert(Signal));

        %Get the spline of the upper envelope
        U=spline(t,Sig);
        env=mkpp(U.breaks,U.coefs);

end

function File_Create(File)

if ~exist(File,'dir')
    mkdir(File);
end

end

function Print(Tag,ode,ode_env,t,File,i,j)

figure()

plot(t,ode.*100,t,ode_env.*100);
title(Tag+" muscle activation from ODE solver")
legend("Windowing pipeline","New pipeline")
xlabel("Time (s)")
ylabel("Muscle activation (%)")
saveas(gcf,File+Tag+" Muscle_activation_ODE_Par_"+i+"_Cyl_"+j+".png");

end

function PrintR(Tag,ode,ode_env,t,File,i,j)

figure()

plot(t,ode_env./ode);
title("Ratio of "+Tag+" muscle activation from ODE solver")
xlabel("Time (s)")
ylabel("Muscle activation Ratio")
saveas(gcf,File+Tag+" Muscle_activation_ODE_Ratio_Par_"+i+"_Cyl_"+j+".png");

end


function [M,sdM,SM,sdSM,all_diffs]=intra_Boxplot(act,sact,i,Tag,File)

if i==1
    J=1:5;
elseif i==2
    J=6:10;
elseif i==3
    J=11:15;
elseif i==4
    J=15:20;
elseif i==5
    J=21:24;
elseif i==6
    J=26:30;
elseif i==7
    J=31:35;
elseif i==8
    J=36:40;
elseif i==9
    J=41:45;
else
    J=46:50;
end

if i~=5
a=act{J(1)};
b=act{J(2)};
c=act{J(3)};
d=act{J(4)};
e=act{J(5)};

sa=sact{J(1)};
sb=sact{J(2)};
sc=sact{J(3)};
sd=sact{J(4)};
se=sact{J(5)};

else
a=act{J(1)};
b=act{J(2)};
c=act{J(3)};
d=act{J(4)};
e=NaN;  

sa=sact{J(1)};
sb=sact{J(2)};
sc=sact{J(3)};
sd=sact{J(4)};
se=NaN;
end

n=max([length(a),length(b),length(c),length(d),length(e)]);

if length(a)<n
    a(end+1:n)=NaN;
    sa(end+1:n)=NaN;
end
if length(b)<n
    b(end+1:n)=NaN;
    sb(end+1:n)=NaN;
end
if length(c)<n
    c(end+1:n)=NaN;
    sc(end+1:n)=NaN;
end
if length(d)<n
    d(end+1:n)=NaN;
    sd(end+1:n)=NaN;
end
if length(e)<n
    e(end+1:n)=NaN;
    se(end+1:n)=NaN;
end

if i~=5
    M=mean([a,b,c,d,e],'omitnan');
    sdM=std([a,b,c,d,e],'omitnan');
    SM=mean([sa,sb,sc,sd,se],'omitnan');
    sdSM=std([a,b,c,d,e],'omitnan');
else
    M=mean([a,b,c,d],'omitnan');
    sdM=std([a,b,c,d],'omitnan');
    SM=mean([sa,sb,sc,sd],'omitnan');
    sdSM=std([a,b,c,d],'omitnan');
end


figure('Position',[100, 100, 900,500])

if i~=5
    data = {[a,sa],[b,sb],[c,sc],[d,sd],[e,se]};
    SS={"Trial 1","Trial 2","Trial 3","Trial 4","Trial 5"};
    colors=["r",'k','b','m', 'g'];
else
    data = {[a,sa],[b,sb],[c,sc],[d,sd]};
    SS={"Trial 1","Trial 2","Trial 3","Trial 4"};
    colors = ["r", "k", "b", "m"]; 
end
all_diffs = []; % before loop
boxplotGroup(data, 'interGroupSpace', 1,'Colors',colors, 'primarylabels', SS,'secondarylabels',{'New pipeline','Existing pipeline'});
ylabel("Muscle activation (%)")
title("Boxplots of "+Tag+" muscle activation of participant "+i);
for t = 1:length(data)
    % For each trial, compare the two groups
    [p_trial, ~] = signrank(data{t}(:,1), data{t}(:,2));
    diffs=data{t}(:,1) - data{t}(:,2);
    all_diffs = [all_diffs, diffs];
    
    if p_trial>1e-3
        pval_str = sprintf('Trial %d p = %.4f', t, p_trial);
    else
        pval_str = sprintf('Trial %d p < %.3f', t, 0.001);
    end
    % Place each text box lower down
    annotation('textbox', [0.15, 0.8-0.05*t, 0.3, 0.05], 'String', pval_str, ...
        'FitBoxToText', 'on', ...
        'EdgeColor', 'none', 'FontSize', 9, 'Color', colors(t));
end

% s=input("You want to save: ","s");
% if s=='y' || s=='Y'
% exportgraphics(gcf,File+"subj_"+i+"_"+Tag+".png")
%end
end


function [b,a,f_cut]=HHT_Filter_Cut_off(hs,ts,fn,fs,HHT_View,ff_HHT,i,j,flag,sig)

%% 1) Convert to power spectrum (Hilbert spectrum)
P = abs(hs).^2;                % power vs time & frequency

%% 2) Marginal spectrum (integrate power over time)
% use trapezoidal integration over time
M = full( trapz(ts, P, 2) );           % column vector, same length as fs

%% 3) Cumulative energy vs frequency
cumE      = cumsum(M);         
cumE_norm = cumE ./ cumE(end); % 0 → 1

%% 4) Choose target fraction of total power (e.g. 95%)
eta = 0.95;                    % keep 95% of power
idx = find(cumE_norm >= eta, 1, 'first');
f_cut = fs(idx);

[b,a]=butter(4,[10 f_cut]./fn,'bandpass');

if isempty(HHT_View)~=1

figure()
hht(emd(sig),2*fn)
hold on;
yline(f_cut,'--','Upper Cut off frequency');
yline(10,'--','Lower Cut off frequency');
hold off;
saveas(gcf,ff_HHT+flag+"_HHT_Par_"+i+"_Cyl_"+j+".png");


end

end


