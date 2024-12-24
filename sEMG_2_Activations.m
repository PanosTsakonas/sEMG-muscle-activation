clc;
close all;
clear all;

fsEMG=1000;
fn=fsEMG/2;

file=getenv('USERPROFILE')+"\OneDrive - University of Warwick\PhD\Hand Trials\Results\Cylindrical Grasp\sEMG Data\";

%Set 1 for zero initial condition in the ODE solver or to [] for non zero
Zero=1;

if isempty(Zero)==0

    ff=file+"Plots\Activation_zero_init\";
    ffr=file+"Plots\Ratio_zero_init\";
else
    ff=file+"Plots\Activation_nonzero_init\";
    ffr=file+"Plots\Ratio_nonzero_init\";
end


%Set to a number to allow the HHT plot of sEMGs
HHT_View=1;

if isempty(HHT_View)==0
    ff_HHT=file+"Plots\HHT\";
end

%Set the filter

[bemg,aemg]=butter(4,[10 450]./fn,'bandpass');

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


for i=1:10

    for j=1:5
         t=[];
        if (i==5 && j==5)
       
        
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
        EDCf=abs(filtfilt(bemg,aemg,EDC))./MVC(i,1);
        FDPf=abs(filtfilt(bemg,aemg,FDP))./MVC(i,3);
        FDSf=abs(filtfilt(bemg,aemg,FDS))./MVC(i,2);
        
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

        % pedc=spline(t,smooth_EDC);
        % pfdp=spline(t,smooth_FDP);
        % pfds=spline(t,smooth_FDS);
        % 
        % ppedc=mkpp(pedc.breaks,pedc.coefs);
        % ppfdp=mkpp(pfdp.breaks,pfdp.coefs);
        % ppfds=mkpp(pfds.breaks,pfds.coefs);

        %Solve the ODE from Blana et al.

        % [~,a_edc_ode]=ode45(@(t,a) activation(t,a,ppedc),t,0);
        % [~,a_fdp_ode]=ode45(@(t,a) activation(t,a,ppfdp),t,0);
        % [~,a_fds_ode]=ode45(@(t,a) activation(t,a,ppfds),t,0);


        figure()
        plot(t,smooth_EDC,t,ppval(EDC_env,t),'--','LineWidth',1)
        hold on;
        plot(t,EDCf);
        hold off;
        title("EDC signal");
        legend("Existing pipeline","New pipeline","EMG signal");
        xlabel("Time (s)");
        ylabel("Normalised EMG signal");
        saveas(gcf,file+"Plots\Signal\EDC_Signal_Par_"+i+"_Cyl_"+j+".png");

        figure()
        plot(t,smooth_FDS,t,ppval(FDS_env,t),'--','LineWidth',1)
        hold on;
        plot(t,FDSf);
        hold off;
        title("FDS signal");
        legend("Existing pipeline","New pipeline","EMG signal");
        xlabel("Time (s)");
        ylabel("Normalised EMG signal");
        saveas(gcf,file+"Plots\Signal\FDS_Signal_Par_"+i+"_Cyl_"+j+".png");

        figure()
        plot(t,smooth_FDP,t,ppval(FDP_env,t),'--','LineWidth',1)
        hold on;
        plot(t,FDPf);
        hold off;
        title("FDP signal");
        legend("Existing pipeline","New pipeline","EMG signal");
        xlabel("Time (s)");
        ylabel("Normalised EMG signal");
        saveas(gcf,file+"Plots\Signal\FDP_Signal_Par_"+i+"_Cyl_"+j+".png");

        figure()
        plot(t,smooth_EDC,t,a_edc_ode_env)
        title("EDC muscle activation");
        legend("Windowing pipeline","New pipeline");
        xlabel("Time (s)");
        ylabel("Muscle activation");
        saveas(gcf,ff+"EDC_Activation_Par_"+i+"_Cyl_"+j+".png");

        figure()
        plot(t,smooth_FDP,t,a_fdp_ode_env)
        title("FDP muscle activation");
        legend("Windowing pipeline","New pipeline");
        xlabel("Time (s)");
        ylabel("Muscle activation");
        saveas(gcf,ff+"FDP_Activation_Par_"+i+"_Cyl_"+j+".png");

        figure()
        plot(t,smooth_FDS,t,a_fds_ode_env)
        title("FDS muscle activation");
        legend("Windowing pipeline","New pipeline");
        xlabel("Time (s)");
        ylabel("Muscle activation");
        saveas(gcf,ff+"FDS_Activation_Par_"+i+"_Cyl_"+j+".png");

         figure()
        plot(t,a_edc_ode_env./smooth_EDC)
        title("Ratio between new and existing pipelines for EDC muscle");
        xlabel("Time (s)");
        ylabel("Ratio of muscle activations");
        saveas(gcf,ffr+"EDC_Ratio_Par_"+i+"_Cyl_"+j+".png");

        figure()
        plot(t,a_fdp_ode_env./smooth_FDP)
        title("Ratio between new and existing pipelines for FDP muscle");
        xlabel("Time (s)");
        ylabel("Ratio of muscle activations");
        saveas(gcf,ffr+"FDP_Ratio_Par_"+i+"_Cyl_"+j+".png");

        figure()
        plot(t,a_fds_ode_env./smooth_FDS)
        title("Ratio between new and existing pipelines for FDS muscle");
        xlabel("Time (s)");
        ylabel("Ratio of muscle activations");
        saveas(gcf,ffr+"FDS_Ratio_Par_"+i+"_Cyl_"+j+".png");

        end
        clear t EDC FDP FDS

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
