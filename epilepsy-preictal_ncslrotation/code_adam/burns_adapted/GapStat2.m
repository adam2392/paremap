%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Original File written by Samuel Burns
% Adapted by: Adam Li
%
% Ver.: 1.0 - Date: 08/14/2015
%
% Description: To find number of clusters using Gap statistics for
% interictal, preictal, seizure and postictal. And perform custom K-mean to
% cluster data.
% 
% Output:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Clearing vars and setup
clear data d1 d2 NumT NumC SimVecs SimWk NumClus *logSimWk 
clear dU dS dV Dmax Dmin Dtmax Dtmin dataname
clear opt

clear
close all

addpath('/home/adamli/MATLAB/code_adam/burns_adapted')

addpath('/Users/adam2392/Documents/MATLAB/Johns Hopkins/epilepsy/code_adam/burns_adapted/MatFiles');
addpath('/Users/adam2392/Documents/MATLAB/Johns Hopkins/epilepsy/code_adam/burns_adapted');
%% Setting up the script's settings
opt = statset('MaxIter',5000);
Patient = 'PY04N008';
dataname = 'AllSZ';              %array to cluster\
FreqBands = {'alpha', 'beta', 'gamma', 'high'};
timeWindow = '600';

% addpath(strcat('/home/adamli/MATLAB/code_adam/burns_adapted/focus_electrode_analysis/', FreqBand))
addpath(['/home/adamli/MATLAB/code_adam/burns_adapted/EVC'])

%%%%% Change if don't want to rescale
rescaled = 0;
szindex = {'8', '20'};
% szindex = {'61', '74', '81'};
name = 'pre';
if rescaled == 0
   szindex = 1; 
end


for f=1:length(FreqBands)
    FreqBand = FreqBands{f}
    
    for i=1:length(szindex)
        cd('/home/adamli/MATLAB/code_adam/burns_adapted')

        mattype = 'q'; % q or pwr_q
        if strcmp(mattype, 'pwr_q')
            mat = 'pwr_';
        else
            mat = '';
        end


        tic;
        %%%%%% Load in EVC files from /EVC/pwr or /EVC/chr
        % eval(['load q',dataname, '_', Patient]) % original chr_matrix
        disp([mat, 'q', dataname, '_', Patient, '_', FreqBand, timeWindow])
        eval(['load ', mat, 'q', dataname, '_', Patient, '_', FreqBand, timeWindow])
        eval(['data = q',dataname,';'])  %set data = qAllII;

        %%%% Try to Rescale and Compare 
        if rescaled
        %     data = resample(data, 1000, 3 * str2num(timeWindow));
            data = load(['q' szindex{i} name '_' Patient '_' timeWindow]);
            eval(['data = data.q', szindex{i}, name, ';']);
            data = resample(data, 3600, str2num(timeWindow));
        end

        %%%% Try to do SVD normally, else do 'economic' setting
        try
            [dU,dS,dV] = svd(data);          %SVD(x) -> USV = X
        catch error
            disp(error)
        %     data = downsample(data, 150);
            data = data(500:3500,:);

            try
                [dU,dS,dV] = svd(data);   %economic version of svd = dU and dS are smaller
            catch error
                [dU,dS,dV] = svd(data,0);   %economic version of svd = dU and dS are smaller
            end

        %     [dU,dS,dV] = svd(data,0);   %economic version of svd = dU and dS are smaller

        end
        tdata = data*dV;            %**** Ask Pierre***      

        [NumT,NumC] = size(data);        %# rows(NumT), cols(NumC) of data(e.g. qallII)

        % %%%% redefine data to make computation work
        % data = data(1:NumT/6,:);
        % [dU,dS,dV] = svd(data);          %SVD(x) -> USV = X   
        % tdata = data*dV;            %**** Ask Pierre***  

        NumSims = 50;                    %# of simulations to run
        NumClus = 5;                     %largest # of clusters to produce

        % loop through columns of data
        for hh = 1:NumC
          Dmax(hh) = max(data(:,hh));
          Dmin(hh) = min(data(:,hh));

          Dtmax(hh) = max(tdata(:,hh));
          Dtmin(hh) = min(tdata(:,hh));
        end

        disp(['Running through patient ', Patient]);
        disp(['Looking at the data: ', dataname]);

        % pool = parpool;             % invoke workers
        % stream = RandStream('mlfg6331_64'); % random number stream
        % options = statset('UseParallel', 1, 'UseSubStreams', 1, 'Streams', stream);

        %% Conduct 'NumSims' # of simulations - Generate 'B' Monte Carlo reference data sets 
        for jj = 1:NumSims
          disp(strcat('On simulation:', ' ', num2str(jj)));
          clear SimVecs WkAll tSimVecs

          disp(['Looping through # of columns in data: ' num2str(NumC)]);
          for hh = 1:NumC
            tSimVecs(1:NumT,hh) = Dtmin(hh) + rand(1,NumT)*(Dtmax(hh)-Dtmin(hh));
          end

          SimVecs = tSimVecs*dV';
          clear tSimVecs

          disp(['Looping through # of rows (timepoints) in data: ', num2str(NumT)]);
          for hh = 1:NumT
            SimVecs(hh,:) = SimVecs(hh,:)/norm(SimVecs(hh,:));
          end

          % generate # of clusters?
          for kk = 1:NumClus
            disp(['Generating ', num2str(kk), ' clusters']);  
            %%%% perform k means with varying # of clusters w/ max iterations set above
            %kk
            clear IDX Centroids Wk
            %[IDX,Centroids] = kmeans(SimVecs,kk,'Options',opt);
            %[IDX,Centroids] = kmeans(SimVecs,kk,'Distance','cityblock','Options',opt);
            [IDX,Centroids] = kmeans(SimVecs,kk,'Distance','cosine','Options',opt);

            Wk = 0;

            % loop through each cluster
            for ff = 1:kk
              clear FFvecs sumW FFlength
              sumW = 0;
              FFlength = length(find(IDX == ff));   % # of observations of cluster 'ff' n_r
              FFvecs = SimVecs(find(IDX == ff),:);

              % sum of all pairwise distances for all points in cluster 'ff'
              disp(['FFlength is : ', num2str(FFlength)]);
              disp(['FFvecs length is : ', num2str(size(FFvecs))]);

        %       tic;
        %       pairdistmatrix = zeros(FFlength,FFlength);
        %       for bb = 1:FFlength
        %           if mod(bb, 1000) == 0
        %               bb
        %           end
        %           for cc=bb+1:FFlength
        %               pairdistmatrix(bb,cc) = norm(FFvecs(bb,:) - FFvecs(cc,:));
        %               pairdistmatrix(cc,bb) = pairdistmatrix(bb,cc);
        %           end
        %       end
        %       toc;

              for bb = 1:FFlength
                  for cc = 1:FFlength   % go through and add 2-norm (xij-xnm)^2
                    sumW = sumW + norm(FFvecs(bb,:)-FFvecs(cc,:));
                  end
              end

              % calculate the pooled within-cluster sum of squares around cluster means
              Wk = Wk + 1/(2*FFlength)*sumW;    
            end

            WkAll(kk) = Wk; % add to array of all W_k
          end

          SimWk(jj,1:NumClus) = WkAll;
        end

        logSimWk = log(SimWk);                  % compute log of W_k
        mlogSimWk = mean(logSimWk);             % compute avg(log(W_k))
        stdlogSimWk = std(logSimWk);            % compute std(log(W_k))
        SimSk = stdlogSimWk*sqrt(1+1/NumSims);  % s_k = sqrt(1 + 1/B) * std(k) | B = # Monte Carlo simulations

        clear s1 s2 pp SZstart SZstop KmeanStates SZlength
        clear opt Wk WkAll
        opt = statset('MaxIter',5000);
        SZstart = 1;
        SZstop = NumT;
        SZlength = NumT;
        s1 = NumT;            

        for hh = 1:s1
          IP_LeadVecs(hh)  = sum(data(hh,:).*data(5,:));
        end

        for hh = 1:s1-1
          dIP_LeadVecs(hh)   = (IP_LeadVecs(hh+1)  - IP_LeadVecs(hh));
        end

        dIP_LeadVecs(s1)   = dIP_LeadVecs(s1-1);
        dIntLeadVecs2c = dIP_LeadVecs; 


        clear dSTD NumSTDs CAxis*
        clear top5 bottom5 TrimSTD sdIntLeadVecs2c
        sdIntLeadVecs2c = sort(dIntLeadVecs2c);
        bottom5 = round(0.025*SZlength);
        top5    = round(0.975*SZlength);
        TrimSTD = std(sdIntLeadVecs2c(bottom5:top5));
        dSTD2 = TrimSTD;
        NumSTDs = 1;

        clear TotalStates StateRun StateTimes NumSTDs* NC KmeanStates 
        clear Silho* Dall* sumd*

        %% Loop through # of clusters we can have
        for KmeanStates = 1:NumClus
        disp(strcat('KmeanStates is: ', num2str(KmeanStates)));
        NC = KmeanStates;
          TotalStates = -1;
          NumSTDs = 1;
          NumIter = 1;

          while TotalStates ~= KmeanStates
            clear NumStates StateRun StateTimes SZtime TotalStates
            StateTimes(1,1:2) = 0;
            NumStates = 0;
            SZtime    = 1;

            if      abs(dIntLeadVecs2c   (1)) >= NumSTDs*dSTD2
              StateRun  = 0;
            elseif (abs(dIntLeadVecs2c(1)) <  NumSTDs*dSTD2) & ...
                   (abs(dIntLeadVecs2c(2)) <  NumSTDs*dSTD2)
              StateRun  = 1;
              NumStates = 1;  
              StateTimes(1,1) = SZstart;
            elseif (abs(dIntLeadVecs2c(1)) <  NumSTDs*dSTD2) & ...
                   (abs(dIntLeadVecs2c(2)) >= NumSTDs*dSTD2)
              StateRun  = 0;
            end

            for gg = 2:SZlength
              SZtime = SZtime + 1;

              if (StateRun == 1) & abs(dIntLeadVecs2c(gg)) >= NumSTDs*dSTD2
                StateTimes(NumStates,2) = SZtime;
                StateRun  = 0;
              end

              if (StateRun == 0) & abs(dIntLeadVecs2c(gg)) < NumSTDs*dSTD2
                StateRun = 1;
                NumStates = NumStates + 1;
                StateTimes(NumStates,1) = SZtime;
              end

              if (dIntLeadVecs2c(gg)   >= NumSTDs*dSTD2) & ...
                 (dIntLeadVecs2c(gg-1) <= -NumSTDs*dSTD2)
                NumStates = NumStates + 1;
                StateTimes(NumStates,1:2) = SZtime;
              end

              if (dIntLeadVecs2c(gg)   <= -NumSTDs*dSTD2) & ...
                 (dIntLeadVecs2c(gg-1) >= NumSTDs*dSTD2)
                NumStates = NumStates + 1;
                StateTimes(NumStates,1:2) = SZtime;
              end

              if (gg == SZlength) & (StateTimes(NumStates,2) == 0)
                StateTimes(NumStates,2) = SZstop;
              end
            end

            TotalStates = NumStates;

            if TotalStates ~= KmeanStates
              NumSTDs = NumSTDs + (TotalStates - KmeanStates)/(log(NumIter+1)*50);
            end

            NumIter = NumIter + 1;
          end

          clear StateVecs StateRun
          for ff = 1:NumStates
            if (StateTimes(ff,2)- StateTimes(ff,1)) > 0
              StateVecs(ff,1:NumC) = ... 
                mean(data(StateTimes(ff,1):StateTimes(ff,2),1:NumC));
            elseif (StateTimes(ff,2)- StateTimes(ff,1)) == 0
              StateVecs(ff,1:NumC) = ... 
                data(StateTimes(ff,1),1:NumC);
            end
          end


          clear IDX Centroids Wk
          %[IDX,Centroids] = kmeans(data,KmeanStates,'start',StateVecs,'Options',opt);
          %[IDX,Centroids] = kmeans(data,KmeanStates,'start',StateVecs,'Distance','cityblock','Options',opt);
          [IDX,Centroids] = kmeans(data,KmeanStates,'start',StateVecs,'Distance','cosine','Options',opt);
          Wk = 0;

          % loop through all the states
          for ff = 1:KmeanStates
            clear FFvecs sumW FFlength
            sumW = 0; %compute successive L1-norm of EVC vectors
            FFlength = length(find(IDX == ff));
            FFvecs   = data(find(IDX == ff),:);

            % sum of all pairwise distances for all points in cluster 'ff'
            for bb = 1:FFlength
                for cc = 1:FFlength
                  sumW = sumW + norm(FFvecs(bb,:)-FFvecs(cc,:));
                end
            end

            Wk = Wk + 1/(2*FFlength)*sumW;
          end

            WkAll(KmeanStates) = Wk;
        end %end loop through K-means

        clear *DataWk
        DataWk(1:NumClus) = WkAll;
        logDataWk(1:NumClus) = log(DataWk);

        %%% Calculate the 'Gap Statistic' (Step 2) 
        clear Gap dGap Qclusters
        Gap = mlogSimWk - logDataWk;    % Gap(k) = mean(log(simWk)) - log(Wk)

        for bb = 1:NumClus-1
          dGap(bb) = Gap(bb) - Gap(bb+1) + SimSk(bb+1); % generate # for inequality in Step 3 of Gap Stat paper
        end

        % plotting the gap statistic graph
        figure(12)
        hold on
        plot(dGap,'r.-')
        grid on

        % find the optimal # of clusters based on gap statistic (Step 3)
        Qclusters = min(find(dGap(2:end) >= 0)) + 1;
        eval(['m',dataname,' = Qclusters;'])

        %% Settings to conduct 'real' k-means
        KmeanStates = Qclusters;
        NC = KmeanStates;
        TotalStates = -1;
        NumSTDs = 1;
        NumIter = 1;

        disp(strcat('Optimal # of clusters based on Gap Statistic is : ', num2str(KmeanStates)));

        %%%% while, the # states, does not equal optimal # of clusters loop
        %%%% do some fancy settings ...?
        while TotalStates ~= KmeanStates        
            clear NumStates StateRun StateTimes SZtime TotalStates
            StateTimes(1,1:2) = 0;
            NumStates = 0;
            SZtime    = 1;

            if      abs(dIntLeadVecs2c   (1)) >= NumSTDs*dSTD2
              StateRun  = 0;
            elseif (abs(dIntLeadVecs2c(1)) <  NumSTDs*dSTD2) & ...
                   (abs(dIntLeadVecs2c(2)) <  NumSTDs*dSTD2)
              StateRun  = 1;
              NumStates = 1;  
              StateTimes(1,1) = SZstart;
            elseif (abs(dIntLeadVecs2c(1)) <  NumSTDs*dSTD2) & ...
                   (abs(dIntLeadVecs2c(2)) >= NumSTDs*dSTD2)
              StateRun  = 0;
            end

            for gg = 2:SZlength
              SZtime = SZtime + 1;

              if (StateRun == 1) & abs(dIntLeadVecs2c(gg)) >= NumSTDs*dSTD2
                StateTimes(NumStates,2) = SZtime;
                StateRun  = 0;
              end

              if (StateRun == 0) & abs(dIntLeadVecs2c(gg)) < NumSTDs*dSTD2
                StateRun = 1;
                NumStates = NumStates + 1;
                StateTimes(NumStates,1) = SZtime;
              end

              if (dIntLeadVecs2c(gg)   >= NumSTDs*dSTD2) & ...
                 (dIntLeadVecs2c(gg-1) <= -NumSTDs*dSTD2)
                NumStates = NumStates + 1;
                StateTimes(NumStates,1:2) = SZtime;
              end

              if (dIntLeadVecs2c(gg)   <= -NumSTDs*dSTD2) & ...
                 (dIntLeadVecs2c(gg-1) >= NumSTDs*dSTD2)
                NumStates = NumStates + 1;
                StateTimes(NumStates,1:2) = SZtime;
              end

              if (gg == SZlength) & (StateTimes(NumStates,2) == 0)
                StateTimes(NumStates,2) = SZstop;
              end
            end

            TotalStates = NumStates;

            if TotalStates ~= KmeanStates
              NumSTDs = NumSTDs + (TotalStates - KmeanStates)/(log(NumIter+1)*50);
            end

            NumIter = NumIter + 1;
        end


        clear StateVecs StateRun
        for ff = 1:NumStates
            if (StateTimes(ff,2)- StateTimes(ff,1)) > 0
              StateVecs(ff,1:NumC) = ... 
                mean(data(StateTimes(ff,1):StateTimes(ff,2),1:NumC));
            elseif (StateTimes(ff,2)- StateTimes(ff,1)) == 0
              StateVecs(ff,1:NumC) = ... 
                data(StateTimes(ff,1),1:NumC);
            end
        end


        % delete(gcp);

        eval(['clear IDX',dataname,' Cents',dataname])
        %eval(['[IDX',dataname,',Cents',dataname,',sumD',dataname,',D',dataname,'] = kmeans(data,Qclusters,''start'',StateVecs,''Options'',opt);'])
        %eval(['[IDX',dataname,',Cents',dataname,',sumD',dataname,',D',dataname,'] = kmeans(data,Qclusters,''start'',StateVecs,''Distance'',''cityblock'',''Options'',opt);'])
        eval(['[IDX',dataname,',Cents',dataname,',sumD',dataname,',D',dataname,'] = kmeans(data,Qclusters,''start'',StateVecs,''Distance'',''cosine'',''Options'',opt);'])

        for bb = 1:Qclusters
          eval(['Cents',dataname,'(bb,:) = Cents',dataname,'(bb,:)/norm(Cents',dataname,'(bb,:));'])
        end

        toc;

        %% Saving the Results of Clustering Into Sep. Mat File

        cd cluster
        if rescaled
            mat = strcat(szindex{i}, 'rescaled_');
        end
        eval(['save ', mat, 'cluster', dataname, '_', Patient, '_', FreqBand, timeWindow, ' IDX', dataname, ' Cents',dataname])

    end
end