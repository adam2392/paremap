%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script analysisSeizeDistribution.m
%
% Description: Analyze the seizure distribution during seizure dynamics of
% each of the seizure. Create dynamic histogram to show the transition
% times during each seizure.
% 
% Output:   no output returned.
%
% Author: Adam Li
%
% Ver.: 1.0 - Date: 10/07/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear
%% 01: Create Settings 
Patient = 'PY04N008';
FreqBand = 'high';
FreqBands = {'alpha', 'beta', 'gamma', 'high'};
% szindices = [70 73 77 79];
% szindices = [61 74 81]; %PY04N007
% szindices = {'124A', '124B', '125'}; %PY05N004
% szindices = {'762' '767' '777' '780' '783' '801'}; %PY04N013
szindices = {'8', '20'}; %PY04N008 
debug_on = 1;

mattype = 'q'; % q, or pwr_q
if strcmp(mattype, 'pwr_q')
    mat = 'pwr_';
else
    mat = '';
end

% load in dataset
filename = strcat(Patient, '_seizureDistribution');
load(filename);

%% 02: loop through each frequency band and seizure
for i=1:length(szindices)
    seizurenum = szindices{i};
    
    disp(['On seizure event #: ' num2str(i)])
    
    % Create figure
    A = figure(1);
    b = subplot(length(szindices), 1, i);
    xlim([0 250])
    ax = gca;
    ax.YTick = [0 1 2 3 4];
    title(['transition time histograms for seizure event: ' num2str(i)])
    xlabel('seconds')
    ylabel('# of occurrences')
    hold on;
    
    % get all the transitions
    for j=1:length(FreqBands)
        FreqBand = FreqBands{j};
        % load in the corresponding seizure transition times in this band
        eval(['transitiontimes = ' FreqBand '.seizure_' seizurenum '.transitiontimes;'])
        
        % store all the transitions for that seizure in a cell array
        transitions{j} = transitiontimes;
        
        len(j) = length(transitiontimes);
    end
    %% 03: Analyze Transition Cell Array
    longestlen = max(len);
    shortestlen = min(len);
    histdata = [];
    colors = {'b', 'r', 'g', 'c', 'y', 'k', 'm', 'w'};
    
    %% Case 1: the histograms for the state transitions that happen in all the bands
    for l=1:shortestlen
       tempdata = cellfun(@(c) c(l), transitions);
       
       % check if there are outliers?
%        [tempdata, idx, outliers] = deleteoutliers(tempdata, 0.05, 1);
%         if ~isempty(outliers)
%             disp(['outliers on state ' num2str(l) ' were'])
%             disp(outliers)
%             disp(['temp data is: '])
%             disp(tempdata)
%         end
       
        histdata = [histdata; tempdata]; 
       
        % plot histogram
        [n, x] = hist(tempdata);
        h = bar(x,n, colors{l});   
        
    end
    disp('Finished adding in all the state transitions for Case 1.')    
    if debug_on
        histdata
    end
    
    % create a temporary cell array to modify
    temptransitions = transitions;  
    
    %% Case 2: the state transitions afterwards
    %initialize removing which freq. bands
    alpharemoved = 0;
    betaremoved = 0;
    gammaremoved = 0;
    highremoved = 0;
    
    % loop through remaining state transitions
    for l=shortestlen+1:longestlen
        % create a temporary cell array to modify
        temptransitions = transitions;
        
        %% find index of deleted row
        for ll=1:length(temptransitions)
            index(ll) = length(transitions{ll}) < l;
        end
        if index(1) == 1 % alpha removed
            alpharemoved = 1;
            alphastring = ' & with alpha removed';
        end
        if index(2) == 1 % beta removed
            betaremoved = 1;
            betastring = ' & with beta removed';
        end
        if index(3) == 1
            gammaremoved = 1; 
            gammastring = ' & with gamma removed';
        end
        if index(4) == 1
            highremoved = 1;
            highstring = ' & with high removed';
        end
        
        if debug_on
            disp(['We are on the ' num2str(l) ' transition'])
        end
        
        %% Continue analysis
        % remove that freq. band from the transition cell array
        temptransitions = temptransitions(cellfun(@(c) size(c,1) > l-1, temptransitions));
        
        if debug_on
            disp('temptransitions')
            disp(temptransitions)
            
            disp('index: ')
            disp(index)
        end
            
%         try
            % -> then grab the data corresponding to that state transition
            tempdata = cellfun(@(c) c(l), temptransitions);
            
        if debug_on
            disp(tempdata)
        end
            
            % put NaN's where the # should be
            k = 1;              % initialize index that will add data
            
            % get vector of 0's and 1's determining, whether this freq.
            % band was removed or not
            removedindices = [alpharemoved, betaremoved, gammaremoved, highremoved];
            bufferdata = tempdata;      % create a buffer data vector
            
            % loop through each freq band and write either NaN or state
            % transition time
            for indice=1:length(removedindices)
                if removedindices(indice) == 1
                    bufferdata(indice) = NaN;
                else
                    bufferdata(indice) = tempdata(k);
                    k = k+1;
                end
            end
            
            % reset the tempdata
            tempdata = bufferdata;
            clear bufferdata
            histdata = [histdata; tempdata];
%             [tempdata, idx, outliers] = deleteoutliers(tempdata, 0.05, 1);
%             if ~isempty(outliers)
%                 disp(['outliers on state ' num2str(l) ' were: '])
%                 disp(outliers)
%             end
%             
            [n, x] = hist(tempdata);
            h = bar(x,n, colors{l});
%         catch error
%             disp('check case 2')
%             throw(error);
%         end
    end
    % save the histdata
    histogramdata{i} = histdata;
%     eval(['seizure_' num2str(szindices(i)) ' = histdata;'])

    % create legend for the histogram plot
    legstring = {'1st state'; '2nd state'; '3rd state'; '4th state'; '5th state'; '6th state'; '7th state'};
    legend(legstring);
end

filename = strcat(Patient, '_seizureDistribution');
save(filename, FreqBands{:}, 'histogramdata'); 
