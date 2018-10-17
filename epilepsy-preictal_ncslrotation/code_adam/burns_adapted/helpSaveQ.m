function [len, filename] = helpSaveQ(pp, InfoTime, homedir, mattype, FPindex, FreqBand, Patient, NumChans2, mAllAdjMat, mAllAdjSTD, focuselectrodes)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% added mod by Adam Li, to get previous recording if the seizure event
% happens very early, or late, so that the pre, post data is very short
%%% For PY04N007 #74 pre needs the #73 post
    %% *** For looking at focus electrodes/channels only
    % PY04N008:
    if focuselectrodes % if we only want to look at the focused electrodes
        if strcmp(Patient, 'PY04N008')
            chans = [1, 2, 3, 9, 10, 11];
            NumChans2 = 6;
            focuselectrodes = 1;
        end

        % % PY04N013
        if strcmp(Patient, 'PY04N013')
            chans = [5,6,7,8,11,12,13,14,15,16,19,20,21,22,23,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75];
            NumChans2 = 44;
            focuselectrodes = 1;
        end
        %   % PY05N004
        if strcmp(Patient, 'PY05N004')
            Patient
            chans = [80,81,82,83,84,85,88,89,90,91,92,93];
            NumChans2 = 12;
            focuselectrodes = 1;
        end
    end

    %% Regular alg.
    NumWins = 0;
    
    FilePath = char(InfoTime.filename(pp));
    FileDir  = FilePath(1:5);
    eval(['cd ',FileDir])

    hdr          = readhdr('eeg.hdr');
    num_channels = hdr.file_fmt.Numb_chans;
    format_file  = hdr.file_fmt.File_format;
    FileName     = FilePath(6:end);
    NumChans     = num_channels;

    %%% Read in the original part i) adjacency matrix to compute "deviation matrix"
    cd adj_matrices
    clear AdjMat

    if ~strcmp(mattype, 'q')
      disp('not q')
      eval(['AdjMat = read_adj_matrix(''adj_pwr_data_',FilePath(FPindex) ...
           ,'_',FreqBand,'.dat'',',num2str(NumChans),');'])
    else
      eval(['AdjMat = read_adj_matrix(''adj_chr_data_',FilePath(FPindex) ...
       ,'_',FreqBand,'.dat'',',num2str(NumChans),');'])
    end
    cd ..
    [NumSecs,n1,n2] = size(AdjMat);   % get array dimensions

     %%% Get the channels of interest 
    if sum(Patient == 'PY04N007') == 8 
        if     sum(FileDir(1:5) == 'rec01') == 5
          AdjMat2 = AdjMat(1:NumSecs,[1:6 8:21 24:27 29:31 33 35:40 42:51 54:56 58:63 65:79 81:87],[1:6 8:21 24:27 29:31 33 35:40 42:51 54:56 58:63 65:79 81:87]);
        end
    end
    
    if sum(Patient == 'PY04N008') == 8 
        if     sum(FileDir(1:5) == 'rec01') == 5
            AdjMat2 = AdjMat(1:NumSecs,[1:40],[1:40]);
        
            %%%%%% Only used when looking at focus electrodes
            if focuselectrodes == 1
                AdjMat2 = AdjMat(1:NumSecs,chans,chans);
            end 
        end
    end
    
    if sum(Patient == 'PY04N012') == 8 
        if     sum(FileDir(1:4) == 'edf0') == 4
          AdjMat2 = AdjMat(1:NumSecs,[1:34 51:98],[1:34 51:98]);
        elseif sum(FileDir(1:4) == 'edf1') == 4
          AdjMat2 = AdjMat;
        elseif sum(FileDir(1:4) == 'edf2') == 4
          AdjMat2 = AdjMat(1:NumSecs,[1:82],[1:82]);
        end
    end
    
    if strcmp(Patient, 'PY04N013')
        AdjMat2 = AdjMat(1:NumSecs,[1 4:10 12:61 63 65 70 72:89],[1 4:10 12:61 63 65 70 72:89]);

        %%%%%% Only used when looking at focus electrodes
        if focuselectrodes == 1
          AdjMat2 = AdjMat(1:NumSecs,chans,chans);
        end 
    end
    
    if sum(Patient == 'PY05N004') == 8 
        if     sum(FileDir(1:4) == 'edf0') == 4
          AdjMat2 = AdjMat(1:NumSecs,[1:12 14 17:37 39:63 65:82 84:100], ...
                                     [1:12 14 17:37 39:63 65:82 84:100]);
        elseif sum(FileDir(1:4) == 'edf1') == 4
          AdjMat2 = AdjMat(1:NumSecs,[1:12 14 17:37 39:81 83:99], ...
                                     [1:12 14 17:37 39:81 83:99]);
        elseif sum(FileDir(1:4) == 'edf2') == 4
          AdjMat2 = AdjMat(1:NumSecs,[1:12 14 17:37 39:81 83:99], ...
                                     [1:12 14 17:37 39:81 83:99]);
        end

        %%%%%% Only used when looking at focus electrodes
        if focuselectrodes == 1
          AdjMat2 = AdjMat(1:NumSecs,chans,chans);
        end
    end
    
    if sum(Patient == 'PY11N009') == 8
        AdjMat2 = AdjMat(1:NumSecs,[1:22 24:30 32:81 85:90 92:105],[1:22 24:30 32:81 85:90 92:105]);
    end

    clear AdjMat
    AdjMat = AdjMat2;
    clear AdjMat2

    % initialize leading EVC vars
    LeadSingVals(1:NumSecs) = 0;
    LeadSingVecs(1:NumSecs,1:NumChans2) = 0;  %rows = seconds, cols = channels

    % loop through each time point (k) in part iii) EVC_k(i)
    for gg = 1:NumSecs
        clear U1 S1 V1 ADJtemp lambda
        lambda = log(2)/2;
        %%% Network Analysis Part ii)
        ADJtemp = squeeze(AdjMat(gg,1:NumChans2,1:NumChans2));
        ADJtemp = (ADJtemp - mAllAdjMat)./mAllAdjSTD;
        ADJtemp = exp(ADJtemp)./(1 + exp(ADJtemp)); % transform to [0,1] to retain "weighted connections" interpretation

        % set diagonals to 0 (electrodes don't connect to each other?)
        for bb = 1:NumChans2
            ADJtemp(bb,bb) = 0;
        end

        ADJtemp(isnan(ADJtemp)) = 0;    % set NaN to 0's

        %%%% perform singular value decomposition
        [U1,S1,V1] = svd(ADJtemp);
        LeadSingVals(gg) = S1(1,1);                             % leading eigenvalues
        LeadSingVecs(gg,1:NumChans2) = abs(U1(1:NumChans2,1));  % leading eigenvectors EVC_k(i)
    end
    cd(homedir)

    [L1,L2] = size(LeadSingVecs);
    eval(['qtest' num2str(pp) '_' Patient FreqBand ' = LeadSingVecs(:,1:NumChans2);']);
    filename = strcat('qtest', num2str(pp), '_', Patient, FreqBand);
    
    save(filename, strcat('qtest', num2str(pp), '_', Patient, FreqBand))

    eval(['len = length(qtest' num2str(pp) '_' Patient FreqBand ');'])
end

