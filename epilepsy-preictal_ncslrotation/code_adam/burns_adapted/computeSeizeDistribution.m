%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script computeSeizeDistribution.m
%
% Description: Compute the seizure distribution during seizure dynamics of
% each of the seizure. Create dynamic histogram to show the transition
% times during each seizure.
%
% Then plot histogram over time
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
Patient = 'PY05N004';
FreqBand = 'high';
FreqBands = {'alpha', 'beta', 'gamma', 'high'};
% szindices = [70 73 77 79];
% szindices = [61 74 81]; %PY04N007
szindices = {'124A', '124B', '125'}; %PY05N004
% szindices = {'762' '767' '777' '780' '783' '801'}; %PY04N013
% szindices = {'8', '20'}; %PY04N008 

% for the k-means settings
szclusters = 6;
distanceFunc = 'cosine'; %or cityblock, correlation, sqeuclidean
options = statset('MaxIter', 3000);

mattype = 'q'; % q, or pwr_q
if strcmp(mattype, 'pwr_q')
    mat = 'pwr_'
else
    mat = ''
end

%% 02: Add Directory Paths and load related Mat files
addpath('/home/adamli/MATLAB/code_adam/burns_adapted/MatFiles')
% addpath(strcat('/home/adamli/MATLAB/code_adam/burns_adapted/cluster/', FreqBand))
addpath(strcat('/home/adamli/MATLAB/code_adam/burns_adapted/EVC/', FreqBand))

addpath(strcat('/Users/adam2392/Documents/MATLAB/Johns Hopkins/epilepsy/code_adam/burns_adapted/EVC/', FreqBand))
% addpath(strcat('/Users/adam2392/Documents/MATLAB/Johns Hopkins/epilepsy/code_adam/burns_adapted/cluster/', FreqBand))
addpath('/Users/adam2392/Documents/MATLAB/Johns Hopkins/epilepsy/code_adam/burns_adapted/MatFiles')

% Temp: Replot the figures
% loop through FreqBands
for j=1:length(FreqBands)
    FreqBand = FreqBands{j}
    addpath(strcat('/home/adamli/MATLAB/code_adam/burns_adapted/cluster/', FreqBand))
    addpath(strcat('/Users/adam2392/Documents/MATLAB/Johns Hopkins/epilepsy/code_adam/burns_adapted/cluster/', FreqBand))
    for i=1:length(szindices)
        seizurenum = szindices{i};
        
        %%%%% Load in preictal cluster files from GapStat.m
        eval(['load ', Patient, '_seizureDistribution'])
        eval(['data = ' FreqBand]);
        eval(['newszidx = data.seizure_' seizurenum '.originalidx;'])
        eval(['medfilteredidx = data.seizure_' seizurenum '.filteredidx;'])
        
        %% Plotting 
        A = figure(j);
        a = subplot(length(szindices), 2, 2*i-1);
        plot(newszidx)
        title([FreqBand sprintf(' band\ntime scaled seizure states for seizure#: ') seizurenum])
        ylim([0 9]);
        ax = gca;
        ax.YTick = [0 1 2 3 4 5 6 7 8 9];
        xlabel('seconds')
        ylabel('seizure states')

        figure(j)
        b = subplot(length(szindices),2,2*i);
        fig = gcf;
        plot(medfilteredidx)
        title([FreqBand sprintf('band\nmedian filtered states for seizure#: ') seizurenum])
        ylim([0 9]);
        ax = gca;
        ax.YTick = [0 1 2 3 4 5 6 7 8 9];
        xlabel('seconds')
        ylabel('seizure states')
    end
end


%% 03: Load in Cluster data to do ictal analysis
% loop through FreqBands
for j=1:length(FreqBands)
    FreqBand = FreqBands{j}
    addpath(strcat('/home/adamli/MATLAB/code_adam/burns_adapted/cluster/', FreqBand))
    addpath(strcat('/Users/adam2392/Documents/MATLAB/Johns Hopkins/epilepsy/code_adam/burns_adapted/cluster/', FreqBand))
    if j==1
        eval(['load ', Patient, '_seizureDistribution'])
    end
    eval(['data = ' FreqBand ';']);
    for i=1:length(szindices)
        seizurenum = szindices{i};
        
        %%%%% Load in preictal cluster files from GapStat.m
        eval(['newszidx = data.seizure_' seizurenum '.originalidx;'])
        eval(['medfilteredidx = data.seizure_' seizurenum '.filteredidx;'])
        
        %% PY04N007
        % alpha
        if strcmp(Patient, 'PY04N007')
            if strcmp(FreqBand, 'alpha')
                statescell = {{8, 4, 8, 6, 5, 7}, ...
                              {8, 4, 8, 6, 5, 7}, ...
                              {8, 4, 8, 6, 7}};
            % beta
            elseif strcmp(FreqBand, 'beta')
                statescell = {{8, 6, 4, 7, 4, 5}, ...
                              {8, 6, 4, 7, 4, 8, 5}, ...
                              {8, 6, 4, 7, 6, 7, 5}};
            % gamma
            elseif strcmp(FreqBand, 'gamma')
                statescell = {{6, 4, 6, 3, 6}, ...
                              {6, 4, 6, 5, 3, 6, 7}, ...
                              {6, 4, 6, 5, 3, 7}};
            % high
            elseif strcmp(FreqBand, 'high')
                statescell = {{3, 6, 4, 5}, ...
                              {3, 6, 4, 5, 7}, ...
                              {3, 6, 4, 5, 7}};
            end
        elseif strcmp(Patient, 'PY04N008')
            if strcmp(FreqBand, 'alpha')
                statescell = {{1,2,1}, ...
                              {1,2,1}};
            % beta
            elseif strcmp(FreqBand, 'beta')
                statescell = {{4, 1, 2}, ...
                              {4, 1, 5, 3}};
            % gamma
            elseif strcmp(FreqBand, 'gamma')
                statescell = {{1, 4, 2, 1}, ...
                              {1, 4, 2, 3, 4}};
            % high
            elseif strcmp(FreqBand, 'high')
                statescell = {{4, 5, 1, 5, 4}, ...
                              {4, 5, 1, 2, 3}};
            end
        elseif strcmp(Patient, 'PY04N013')
            if strcmp(FreqBand, 'alpha')
                statescell = {{1, 2}, ...
                              {2, 1, 2, 1}, ...
                              {2, 1, 2}, ...
                              {2, 1, 2, 1}, ...
                              {2, 1}, ...
                              {1, 2, 1, 2}};
            % beta
            elseif strcmp(FreqBand, 'beta')
                statescell = {{1, 2, 3, 2}, ...
                              {3, 1, 3, 2}, ...
                              {3, 1, 3, 2}, ...
                              {3, 1, 3, 1}, ...
                              {3, 1, 2, 3}, ...
                              {3, 1, 2}};
            % gamma
            elseif strcmp(FreqBand, 'gamma')
                statescell = {{3, 2, 1}, ...
                              {4, 3, 2}, ...
                              {4, 3, 2}, ...
                              {4, 3, 2}, ...
                              {4, 3, 2}, ...
                              {4, 3, 1}};
            % high
            elseif strcmp(FreqBand, 'high')
                statescell = {{2, 3, 1}, ...
                              {2, 3, 1}, ...
                              {2, 3, 1}, ...
                              {2, 3, 2, 1}, ...
                              {2, 3, 2, 1}, ...
                              {2, 3, 1}};
            end
        elseif strcmp(Patient, 'PY05N004')
            if strcmp(FreqBand, 'alpha')
                statescell = {{2, 5, 4, 1}, ...
                              {2, 3, 4, 1}, ...
                              {2, 5, 1, 2, 3, 4, 1}};
            % beta
            elseif strcmp(FreqBand, 'beta')
                statescell = {{4, 2, 5, 3}, ...
                              {4, 2, 5, 3}, ...
                              {2, 5, 3, 2, 5, 2}};
            % gamma
            elseif strcmp(FreqBand, 'gamma')
                statescell = {{1, 2}, ...
                              {1, 2}, ...
                              {1, 2, 1, 2}};
            % high
            elseif strcmp(FreqBand, 'high')
                statescell = {{1, 2}, ...
                              {1, 2}, ...
                              {1, 2, 1, 2}};
            end
        end
        
        %% 04: Identify list of times when state changesuntitled.bmp
        %%%%%% manually code in state vectors to check
        % set the states to the correct index to check
        states = statescell{i};
        windowedidx = medfilteredidx;
        timeofchange = 1;
        offset = 0;
        transitiontimes = zeros(length(states), 1);
        kk = 1;
        
        % select the figure and corresponding seizure subplot
        B = figure(j);
        b = subplot(length(szindices),2,2*i);
        
        while kk ~= length(states)+1
            % first state
            if kk == 1
                timefirststate = find(windowedidx == states{kk}, 1);
                % Case 1: starts on correct state
                if timefirststate == 1
                    disp(['seizure: ' seizurenum ' and ' FreqBand ' is fine.'])
                % Case 2: does not start on correct state, we wanted
                else
                    disp(['seizure: ' seizurenum ' and ' FreqBand ' is not fine.'])
                    offset = timefirststate;                  % store the offset
                    windowedidx = medfilteredidx(offset:end); % create new window
                    
                    % set timefirststate to 1 for handling the offset
                    timefirststate = 1;
                end
                transitiontime = timefirststate;
            % every other state
            else
                % find the transition time of the next state
                % kk=2,3,...length(states)
                transitiontime = find(windowedidx == states{kk}, 1); 
                windowedidx = windowedidx(transitiontime:end); % create new window
                
                % set variables for next loop
                transitiontime = transitiontime + offset -1;
                offset = transitiontime;
                
                % draw line
                line([transitiontime transitiontime], [0 9], 'Color', 'k')
            end
            
            % add time to the transition time vector
            transitiontimes(kk) = transitiontime;
            kk = kk+1;
        end
        % save figure as .fig
        saveas(B, strcat(Patient, '_', FreqBand, '.fig'))
    
        % save the transition times into the corresponding struct inside
        % the (Patient)_seizureDistribution mat file
        eval([FreqBand '.seizure_' seizurenum '.transitiontimes = transitiontimes;']);
    end
end

filename = strcat(Patient, '_seizureDistribution');
save(filename, FreqBands{:});

%% PY04N007
% alpha
% statescell = {{7, 4, 7, 6, 5, 7}, ...
%               {7, 4, 7, 6, 5, 7}, ...
%               {7, 4, 7, 6, 7}};
% beta
% statescell = {{7, 6, 4, 7, 4, 5}, ...
%               {7, 6, 4, 7, 4, 7, 5}, ...
%               {7, 6, 4, 7, 6, 7, 5}};
