clc;
close all;
clear all;

%Set your own sampling frequency
fsEMG=1000;
fn=fsEMG/2;

file=getenv('USERPROFILE')+"Set file path location of sEMG data";

%Set 1 for zero initial condition in the ODE solver or to [] for non zero
Zero=1;

if isempty(Zero)==0

    ff=file+"Plots\Activation_zero_init\";
   
    ffr=file+"Plots\Ratio_zero_init\";
    
else
    ff=file+"Plots\Activation_nonzero_init\";
    ffr=file+"Plots\Ratio_nonzero_init\";
end

%This checks whether the folders exsit that will have all the plots
File_Create(ff);
File_Create(ffr);

ff_Sig=file+"Plots\Signal\";
File_Create(ff_Sig);

ff_ODE=file+"Plots\Activation ODE zero init both\";
ff_ODE_R=file+"Plots\Ratio of Activation ODE zero init both\";
f_box=file+"Plots\BoxPlots\intra_subject\";
f_inter=file+"Plots\BoxPlots\inter_subject\";
File_Create(f_box);
File_Create(ff_ODE);
File_Create(ff_ODE_R);
File_Create(f_inter);

%Set to a number to allow the HHT plot of sEMGs
HHT_View=[];

if isempty(HHT_View)==0
    ff_HHT=file+"Plots\HHT\";
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

        if isempty(HHT_View)~=1

            figure()
            hht(emd(EDC),fsEMG);
            title("HHT of EDC signal");
            saveas(gcf,ff_HHT+"EDC_HHT_Par_"+i+"_Cyl_"+j+".png");
            figure()
            hht(emd(FDS),fsEMG);
            title("HHT of FDS signal");
            saveas(gcf,ff_HHT+"FDS_HHT_Par_"+i+"_Cyl_"+j+".png");
            figure()
            hht(emd(FDP),fsEMG);
            title("HHT of FDP signal");
            saveas(gcf,ff_HHT+"FDP_HHT_Par_"+i+"_Cyl_"+j+".png");
        end
    
        for m=1:length(FDP)

            if m==1
                t(m)=0;
            else
                t(m)=t(m-1)+1/fsEMG;
            end
        end

        %Band pass [10-450], rectify and normalise

        [bEDC,aEDC]=Filter_Cut_Off(EDC,fn);
        [bFDP,aFDP]=Filter_Cut_Off(FDP,fn);
        [bFDS,aFDS]=Filter_Cut_Off(FDS,fn);

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
for i=1:10
[M_EDC,sd_EDC,SM_EDC,sd_SEDC]=intra_Boxplot(edc,sedc,i,"EDC",f_box);
[M_FDS,sd_FDS,SM_FDS,sd_SFDS]=intra_Boxplot(fds,sfds,i,"FDS",f_box);
[M_FDP,sd_FDP,SM_FDP,sd_SFDP]=intra_Boxplot(fdp,sfdp,i,"FDP",f_box);

medc{i}=M_EDC;
sdedc{i}=sd_EDC;
sdsedc{i}=sd_SEDC;
smedc{i}=SM_EDC;
mfds{i}=M_FDS;
sdfds{i}=sd_FDS;
sdsfds{i}=sd_SFDS;
smfds{i}=SM_FDS;
mfdp{i}=M_FDP;
sdfdp{i}=sd_FDP;
sdsfdp{i}=sd_SFDP;
smfdp{i}=SM_FDP;
end


inter_Boxplot(medc,smedc,sdedc,sdsedc,"EDC",f_inter);
inter_Boxplot(mfds,smfds,sdfds,sdsfds,"FDS",f_inter);
inter_Boxplot(mfdp,smfdp,sdfdp,sdsfdp,"FDP",f_inter);




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

function [b,a]=Filter_Cut_Off(sig,fn)

[~,~,~,ifs]=hht(emd(sig),2*fn);

up=max(max(abs(ifs(:,1:4))));


[b,a]=butter(4,[10 up]./fn,'bandpass');


end

function [M,sdM,SM,sdSM]=intra_Boxplot(act,sact,i,Tag,File)

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


% figure('Position',[100, 100, 900,500])
% 
% if i~=5
%     data = {[a,sa],[b,sb],[c,sc],[d,sd],[e,se]};
%     SS={"Trial 1","Trial 2","Trial 3","Trial 4","Trial 5"};
%     colors=["r",'k','b','m', 'g'];
% else
%     data = {[a,sa],[b,sb],[c,sc],[d,sd]};
%     SS={"Trial 1","Trial 2","Trial 3","Trial 4"};
%     colors = ["r", "k", "b", "m"]; 
% end
% 
% boxplotGroup(data, 'interGroupSpace', 1,'Colors',colors, 'primarylabels', SS,'secondarylabels',{'New pipeline','Existing pipeline'});
% ylabel("Muscle activation (%)")
% title("Boxplots of "+Tag+" muscle activation of participant "+i);
% for t = 1:length(data)
%     % For each trial, compare the two groups
%     [p_trial, ~] = signrank(data{t}(:,1), data{t}(:,2));
%     if p_trial>1e-3
%         pval_str = sprintf('Trial %d p = %.4f', t, p_trial);
%     else
%         pval_str = sprintf('Trial %d p < %.3f', t, 0.001);
%     end
%     % Place each text box lower down
%     annotation('textbox', [0.15, 0.8-0.05*t, 0.3, 0.05], 'String', pval_str, ...
%         'FitBoxToText', 'on', ...
%         'EdgeColor', 'none', 'FontSize', 9, 'Color', colors(t));
% end
% s=input("You want to save: ","s");
% if s=='y' || s=='Y'
% exportgraphics(gcf,File+"subj_"+i+"_"+Tag+".png")
% end
end


function inter_Boxplot(M,SM,sdM,sdSM,Tag,File)
F_Fit=File+"fit_"+Tag+"\";
File_Create(F_Fit);
F_Box=File+Tag+" boxplot\";
File_Create(F_Box);
% for i=1:10
% figure()
% 
% errorbar(M{i},SM{i},sdM{i},sdM{i},sdSM{i},sdSM{i},'o');
% ylabel("Muscle activation existing pipeline (%)");
% xlabel("Muscle activation new pipeline (%)");
% hold on;
% P=polyfitZero(M{i},SM{i},1);
% x=[min([M{i}-sdM{i},SM{i}-sdSM{i}]),max([M{i}+sdM{i},SM{i}+sdSM{i}])];
% xn=linspace(x(1),x(2),100);
% plot(xn,polyval(P,xn));
% hold off;
% legend("Data","Y= "+P(1)+" x",'Location','best');
% rs=R2(SM{i},polyval(P,M{i}));
% rms=rmse(SM{i},polyval(P,M{i}));
% title("Errorbars and linear fit of "+Tag+" muscle with R^2: "+rs+" and RMSE: "+rms);
% exportgraphics(gcf,F_Fit+"P_"+i+".png");
% end

M{5}=[M{5},NaN];
SM{5}=[SM{5},NaN];

data={[M{1}',SM{1}'],[M{2}',SM{2}'],[M{3}',SM{3}'],...
[M{4}',SM{4}'],[M{5}',SM{5}'],[M{6}',SM{6}'],...
[M{7}',SM{7}'],[M{8}',SM{8}'],[M{9}',SM{9}'],[M{10}',SM{10}']};

SS={"P1","P2","P3","P4","P5","P6","P7","P8","P9","P10"};

colors = [
    0.1216, 0.4667, 0.7059;   % blue
    1.0000, 0.4980, 0.0549;   % orange
    0.1725, 0.6275, 0.1725;   % green
    0.8392, 0.1529, 0.1569;   % red
    0.5804, 0.4039, 0.7412;   % purple
    0.5490, 0.3373, 0.2941;   % brown
    0.8902, 0.4667, 0.7608;   % pink
    0.4980, 0.4980, 0.4980;   % gray
    0.7373, 0.7412, 0.1333;   % olive
    0.0902, 0.7451, 0.8118    % cyan
];

figure('Position',[100, 100, 900,500])
boxplotGroup(data, 'interGroupSpace', 1,'Colors',colors, 'primarylabels', SS,'secondarylabels',{'New pipeline','Existing pipeline'});
ylabel("mean muscle activation (%)")
title("Boxplots of "+Tag+" muscle activation across participants");
p_wilcoxon = NaN(10,1);
HigherMethod= "";
for p = 1:length(data)
    % Extract data for participant p
    new_pipeline = data{p}(:,1);
    existing_pipeline = data{p}(:,2);

    % Wilcoxon signed-rank test
    p_wilcoxon(p) = signrank(new_pipeline, existing_pipeline);

    if p_wilcoxon(p) < 0.05
        
        if median(new_pipeline, 'omitnan') > median(existing_pipeline, 'omitnan')
            HigherMethod(p) = "New pipeline";
        else
            HigherMethod(p) = "Existing pipeline";
        end
    else
        
        HigherMethod(p) = "-";
    end



    % Format p-value string
    if p_wilcoxon(p) < 0.001
        pval_str = 'p < 0.001';
    else
        pval_str = sprintf('p = %.3f', p_wilcoxon(p));
    end

    % Determine x and y position for annotation
    x_pos = (p-1)*2 + 1.5; % Centered between the two boxes for participant p
    y_pos = max([new_pipeline; existing_pipeline]) * 1.05; % Slightly above the top whisker

    % Add annotation (subtle style)
    text(x_pos, y_pos, pval_str, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'FontSize', 9, ...
        'Color', colors(p,:), ...
        'FontAngle', 'italic', ...
        'Interpreter', 'none');
end
%exportgraphics(gcf,F_Box+"Boxplot_"+Tag+".png");
T_pvals = table((1:10)', p_wilcoxon,HigherMethod', 'VariableNames', {'Participant', 'Wilcoxon_p_value','Higher Activation'});
writetable(T_pvals, F_Box+"Boxplot_"+Tag+"Wilcoxon_pvalues_per_subject.csv");
% Example dimensions (replace with your actual numbers)
N_subjects = 10;
N_trials = 5;
N_methods = 2; % e.g., New and Existing

% Create vectors for each variable
Subject = repelem((1:N_subjects)', N_trials * N_methods);
Trial = repmat(repelem((1:N_trials)', N_methods), N_subjects, 1);
Method = repmat({'New pipeline'; 'Existing pipeline'}, N_subjects * N_trials, 1);
Activation=[];
for i=1:10
    a=M{i};
    b=SM{i};

    for j=1:5

        if i==5 && j==5
            Activation=[Activation;NaN;NaN];
        else
        Activation=[Activation;a(j);b(j)];
        end
    end
end

% Create table
activationTable = table(Subject, Trial, Method, Activation);

lme = fitlme(activationTable, 'Activation ~ Method + Trial + (1|Subject)');
% Export summary tables
% Fixed effects:
fixedEffectsTable = dataset2table(lme.Coefficients);
writetable((fixedEffectsTable), F_Box+"Boxplot_"+Tag+"lme_fixed_effects.csv");

% % Random effects (inter-subject variability):
 randomEffectsTable = randomEffects(lme);
 writetable((randomEffectsTable), F_Box+"Boxplot_"+Tag+"lme_random_effects.csv");

% ANOVA table for model:
anovaTable = dataset2table(anova(lme));
writetable((anovaTable), F_Box+"Boxplot_"+Tag+"lme_anova.csv");


end
