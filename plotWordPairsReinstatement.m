close all
clc
clear all

%%- Section for Pulling reinstatement matrices already produced and just
%%make different plots
ANALYSIS_TYPE = {'within_blocks', 'across_blocks'};
EVENT_SYNC = {'probeon', 'vocalization', 'matchword'};

subj = 'NIH034';
ANALYSIS = ANALYSIS_TYPE{2};
SYNC = EVENT_SYNC{2};

disp(ANALYSIS)
disp(SYNC)

% file dir for all the saved mat files
fileDir = strcat('./Figures/', subj, '/reinstatement_mat/', ANALYSIS, '_', SYNC, '/');

% load in an example data directory to get session names and block number
dataDir = strcat('./condensed_data_', subj);
dataDir = fullfile(dataDir, 'morlet_spec');
sessions = dir(dataDir);
sessions = {sessions(3:end).name};
if strcmp(subj, 'NIH039')
    sessions = sessions([1,2,4]);
elseif strcmp(subj, 'NIH034')
    sessions = sessions([3, 4]);
end
blocks = dir(fullfile(dataDir, sessions{1}));
blocks = {blocks(3:end).name};

% set which blocks to analyze
if strcmp(ANALYSIS, 'across_blocks')
    lenBlocks = length(blocks)-1;
else
    lenBlocks = length(blocks);
end

sessions

avgeReinstatementMat = [];

if strcmp(SYNC, 'vocalization')
    ticks = [6:10:56];
    labels = [-3:1:2];
    timeZero = 36;
    
    labels = [-1.5:0.5:-0.5];
    ticks = [1:5:11];
    
elseif strcmp(SYNC, 'matchword')
    ticks = [6:10:56];
    labels = [-4:1:1];
    timeZero = 46;
elseif strcmp(SYNC, 'probeon')
    ticks = [6:10:56];
    labels = [0:1:5];
    timeZero = 6;
end
% set linethickness
LT = 1.5;

%- LOOP THROUGH SESSIONS
for iSesh=1:length(sessions),
    
    %%- Initialize variables to store averaged reinstatement per session
    avgeSameSessionReinstatement = [];
    avgeDiffSessionReinstatement = [];
    avgeProbeSessionReinstatement = [];
    avgeTargetSessionReinstatement = [];
    avgeReverseSessionReinstatement = [];
    
    % Looking at the entire reinstatement
    % figureFile = strcat('./Figures/', subj, '/reinstatement/', ANALYSIS, '_', SYNC, '/', ...
    %     'summaryttest_', sessions{iSesh});

    % looking at the period before vocalization
    figureDir = strcat('./Figures/', subj, '/preVocalization/', ANALYSIS, '/');
    figureFile = strcat('./Figures/', subj, '/preVocalization/', ANALYSIS, '/', ...
        'summary_', sessions{iSesh});
    %%- LOOP THROUGH BLOCKS
    for iBlock=1:lenBlocks,
        if strcmp(ANALYSIS, 'across_blocks')
            sessionBlockName = strcat(sessions{iSesh}, '-', blocks{iBlock}, 'vs', blocks{iBlock+1});
        else
            sessionBlockName = strcat(sessions{iSesh}, '-', blocks{iBlock});
        end
        fileToLoad = fullfile(fileDir, sessionBlockName);
        data = load(fileToLoad);
        
        %%- LOAD IN THE DATA SAVED 
        eventSame = data.eventSame;
        eventDiff = data.eventDiff;
        featureSame = data.featureSame;
        featureDiff = data.featureDiff;
        if strcmp(ANALYSIS, 'across_blocks')
            eventProbe = data.eventProbe;
            eventReverse = data.eventReverse;
            eventTarget = data.eventTarget;
            featureReverse = data.featureReverse;
            featureTarget = data.featureTarget;
            featureProbe = data.featureProbe;
            
            %%- Section to make averaged responses
            if isempty(avgeProbeSessionReinstatement)
                avgeProbeSessionReinstatement = eventProbe; %permute(eventSame, [3, 1, 2]);
                avgeTargetSessionReinstatement = eventTarget;
                avgeReverseSessionReinstatement = eventReverse;
            else
                avgeProbeSessionReinstatement = cat(1, avgeProbeSessionReinstatement, eventProbe);
                avgeTargetSessionReinstatement = cat(1, avgeTargetSessionReinstatement, eventTarget);
                avgeReverseSessionReinstatement = cat(1, avgeReverseSessionReinstatement, eventReverse);
            end
        end
       
        %%- Section to make averaged responses
        if isempty(avgeSameSessionReinstatement)
            avgeSameSessionReinstatement = eventSame; %permute(eventSame, [3, 1, 2]);
            avgeDiffSessionReinstatement = eventDiff;
        else
            avgeSameSessionReinstatement = cat(1, avgeSameSessionReinstatement, eventSame);
            avgeDiffSessionReinstatement = cat(1, avgeDiffSessionReinstatement, eventDiff);
        end
        
        % rand sample down the different word pair feature mat -> match
        % size
        randIndices = randsample(size(eventDiff,1), size(eventSame,1));
        eventDiff = eventDiff(randIndices,:,:);
        %% Plotting
%         figure
%         subplot(311)
%         imagesc(squeeze    figureFile = strcat('./Figures/', subj, '/reinstatement/', ANALYSIS, '_', SYNC, '/', ...
%         'summaryttest_', sessions{iSesh});(mean(eventSame(:, :, :),1)));
%         title(['Same Pairs Cosine Similarity for Block ', num2str(iBlock-1) ...
%             ' with ', num2str(size(eventDiff,1)), ' events'])
%         hold on
%         xlabel('Time (seconds)');
%         ylabel('Time (seconds)');
%         ax = gca;
%         axis square
%         ax.YTick = ticks;
%         ax.YTickLabel = labels;
%         ax.XTick = ticks;
%         ax.XTickLabel = labels;
%         colormap('jet');
%         set(gca,'tickdir','out','YDir','normal');
%         set(gca, 'box', 'off');
%         colorbar();
%         clim = get(gca, 'clim');
%         hold on
%         plot(get(gca, 'xlim'), [timeZero timeZero], 'k', 'LineWidth', LT)
%         plot([timeZero timeZero], get(gca, 'ylim'), 'k', 'LineWidth', LT)
%        
%         subplot(312);
%         imagesc(squeeze(mean(eventDiff(:, :, :),1)));
%         title(['Different Word Pairs Cosine Similarity for Block ', num2str(iBlock-1)])
%         hold on
%         xlabel('Time (seconds)');
%         ylabel('Time (seconds)');
%         ax = gca;
%         axis square
%         ax.YTick = ticks;
%         ax.YTickLabel = labels;
%         ax.XTick = ticks;
%         ax.XTickLabel = labels;
%         colormap('jet');
%         set(gca,'tickdir','out','YDir','normal');
%         set(gca, 'box', 'off');
%         colorbar();
%         set(gca, 'clim', clim);
%         hold on
%         plot(get(gca, 'xlim'), [timeZero timeZero], 'k', 'LineWidth', LT)
%         plot([timeZero timeZero], get(gca, 'ylim'), 'k', 'LineWidth', LT)
%         
%         subplot(313);
%         imagesc(squeeze(mean(eventSame(:, :, :),1)) - squeeze(mean(eventDiff(:, :, :),1)));
%         title(['Same-Different Word Pairs Cosine Similarity for Block ', num2str(iBlock-1)])
%         hold on
%         xlabel('Time (seconds)');
%         ylabel('Time (seconds)');
%         ax = gca;
%         axis square
%         ax.YTick = ticks;
%         ax.YTickLabel =    figureFile = strcat('./Figures/', subj, '/reinstatement/', ANALYSIS, '_', SYNC, '/', ...
%         'summaryttest_', sessions{iSesh}); labels;
%         ax.XTick = ticks;
%         ax.XTickLabel = labels;
%         colormap('jet');
%         set(gca,'tickdir','out','YDir','normal');
%         set(gca, 'box', 'off');
%         colorbar();
%         hold on
%         plot(get(gca, 'xlim'), [timeZero timeZero], 'k', 'LineWidth', LT)
%         plot([timeZero timeZero], get(gca, 'ylim'), 'k', 'LineWidth', LT)
        
        %%- Save Image
%         print(figureFile, '-dpng', '-r0')
%         savefig(figureFile)
        
        pause(0.1);
    end % loop through blocks
    
    size(avgeSameSessionReinstatement)
    size(avgeDiffSessionReinstatement)
    
    eventSame = avgeSameSessionReinstatement;
    eventDiff = avgeDiffSessionReinstatement;
    
    %%- make a bar plot for the mean/sem at timeZero, timeZero
    if strcmp(ANALYSIS, 'across_blocks')
        eventProbe = avgeProbeSessionReinstatement;
        eventTarget = avgeTargetSessionReinstatement;
        eventReverse = avgeReverseSessionReinstatement;
        
%         eventSame = eventSame(:, timeZero, timeZero);
%         eventDiff = eventDiff(:, timeZero, timeZero);
        
        %%- Make bar plots from 0,0
        sameMean = mean(eventSame(:,timeZero,timeZero), 1);
        sameSem = std(eventSame(:,timeZero,timeZero))/sqrt(size(eventSame,1));
        diffMean = mean(eventDiff(:,timeZero,timeZero), 1);
        diffSem = std(eventDiff(:,timeZero,timeZero))/sqrt(size(eventDiff,1));
        probeMean = mean(eventProbe(:,timeZero,timeZero), 1);
        probeSem = std(eventProbe(:,timeZero,timeZero))/sqrt(size(eventProbe,1));
        targetMean = mean(eventTarget(:,timeZero,timeZero), 1);
        targetSem = std(eventTarget(:,timeZero,timeZero))/sqrt(size(eventTarget,1));
        reverseMean = mean(eventReverse(:,timeZero,timeZero), 1);
        reverseSem = std(eventReverse(:,timeZero,timeZero))/sqrt(size(eventReverse,1));
        
        meanVec = [sameMean, diffMean, reverseMean, probeMean, targetMean];
        semVec = [sameSem, diffSem, reverseSem, probeSem, targetSem];
        figure
        hold on
        p1 = bar(1:5, meanVec);
        title(['Cosine Similarity Across Blocks Reinstatement For Different Word Pair Groups']);
        set(p1, 'FaceColor', 'black');
        p2 = errorbar(1:5, meanVec, semVec, '.');
        set(p2, 'Color', 'red');
        xlabel('Word Pair Groups');
        ylabel('Cosine Similarity at 0,0');
        ax = gca;
        ax.XTick = 1:5;
        ax.XTickLabel = {'Same Pairs', 'Diff Pairs', 'Reverse', 'Probe Overlap', 'Target Overlap'};

        % save the figure
        barFigFile = strcat('./Figures/', subj, '/', sessions{iSesh}, '_summaryBarPlot');
        print(barFigFile, '-dpng', '-r0');
    end

    % change the windows of analysis to before vocalization
    eventSame = eventSame(:, timeZero-10:timeZero, timeZero-10:timeZero);
    eventDiff = eventDiff(:, timeZero-10:timeZero, timeZero-10:timeZero);
    
    figure
    subplot(311)
    imagesc(squeeze(mean(eventSame(:, :, :),1)));
    title(['Same Pairs Cosine Similarity for ', sessions{iSesh}, ...
        ' with ', num2str(size(eventDiff,1)), ' events'])
    hold on
    xlabel('Time (seconds)');
    ylabel('Time (seconds)');
    ax = gca;
    axis square
    ax.YTick = ticks;    
    ax.YTickLabel = labels;
    ax.XTick = ticks;
    ax.XTickLabel = labels;
    colormap('jet');
    set(gca,'tickdir','out','YDir','normal');
    set(gca, 'box', 'off');
    colorbar();
    clim = get(gca, 'clim');
    hold on
%     plot(get(gca, 'xlim'), [timeZero timeZero], 'k', 'LineWidth', LT)
%     plot([timeZero timeZero], get(gca, 'ylim'), 'k', 'LineWidth', LT)

    subplot(312);
    imagesc(squeeze(mean(eventDiff(:, :, :),1)));
    title(['Different Word Pairs Cosine Similarity'])
    hold on
    xlabel('Time (seconds)');
    ylabel('Time (seconds)');
    ax = gca;
    axis square
    ax.YTick = ticks;
    ax.YTickLabel = labels;
    ax.XTick = ticks;
    ax.XTickLabel = labels;
    colormap('jet');
    set(gca,'tickdir','out','YDir','normal');
    set(gca, 'box', 'off');
    set(gca, 'clim', clim);
    colorbar();
    hold on
%     plot(get(gca, 'xlim'), [timeZero timeZero], 'k', 'LineWidth', LT)
%     plot([timeZero timeZero], get(gca, 'ylim'), 'k', 'LineWidth', LT)

    subplot(313);
    imagesc((squeeze(mean(eventSame(:, :, :),1)) - squeeze(mean(eventDiff(:, :, :),1))));
    title(['Same-Different Word Pairs Cosine Similarity'])
    hold on
    xlabel('Time (seconds)');
    ylabel('Time (seconds)');
    ax = gca;
    axis square
    ax.YTick = ticks;
    ax.YTickLabel = labels;
    ax.XTick = ticks;
    ax.XTickLabel = labels;
    colormap('jet');
    set(gca,'tickdir','out','YDir','normal');
    set(gca, 'box', 'off');
%     set(gca, 'clim', [0, 0.05]);
    colorbar();
    hold on
%     plot(get(gca, 'xlim'), [timeZero timeZero], 'k', 'LineWidth', LT)
%     plot([timeZero timeZero], get(gca, 'ylim'), 'k', 'LineWidth', LT)
    
    figureFile
    if ~exist(figureDir)
        mkdir(figureDir);
    end
    print(figureFile, '-dpng', '-r0')
end % loop through sessions