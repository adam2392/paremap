%%- Section for Pulling reinstatement matrices already produced and just
%%make different plots
ANALYSIS_TYPE = {'within_blocks', 'across_blocks'};
EVENT_SYNC = {'probeon', 'vocalization', 'matchword'};

subj = 'NIH034';
ANALYSIS = ANALYSIS_TYPE{1};
SYNC = EVENT_SYNC{2};

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
%- LOOP THROUGH SESSIONS
for iSesh=1:length(sessions),
    %%- LOOP THROUGH BLOCKS
    for iBlock=1:lenBlocks,
        if strcmp(ANALYSIS, 'across_blocks')
            sessionBlockName = strcat(sessions{iSesh}, '-', blocks{iBlock}, 'vs', blocks{iBlock+1});
        else
            sessionBlockName = strcat(sessions{iSesh}, '-', blocks{iBlock});
        end
        fileToLoad = fullfile(fileDir, sessionBlockName);
        data = load(fileToLoad);
        
        eventSame = data.eventSame;
        eventDiff = data.eventDiff;
        featureSame = data.featureSame;
        featureDiff = data.featureDiff;
        
        size(eventSame)
        size(eventDiff)
        
        % rand sample down the different word pair feature mat -> match
        % size
        randIndices = randsample(size(eventDiff,1), size(eventSame,1));
        eventDiff = eventDiff(randIndices,:,:);
        
        if strcmp(SYNC, 'vocalization')
            ticks = [6:10:56];
            labels = [-3:1:2];
            timeZero = 36;
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
        
        %%- Plotting
        figure
        subplot(311)
        imagesc(squeeze(mean(eventSame(:, :, :),1)));
        title(['Same Pairs Cosine Similarity for Block ', num2str(iBlock-1) ...
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
        plot(get(gca, 'xlim'), [timeZero timeZero], 'k', 'LineWidth', LT)
        plot([timeZero timeZero], get(gca, 'ylim'), 'k', 'LineWidth', LT)
       
        subplot(312);
        imagesc(squeeze(mean(eventDiff(:, :, :),1)));
        title(['Different Word Pairs Cosine Similarity for Block ', num2str(iBlock-1)])
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
        set(gca, 'clim', clim);
        hold on
        plot(get(gca, 'xlim'), [timeZero timeZero], 'k', 'LineWidth', LT)
        plot([timeZero timeZero], get(gca, 'ylim'), 'k', 'LineWidth', LT)
        
        subplot(313);
        imagesc(squeeze(mean(eventSame(:, :, :),1)) - squeeze(mean(eventDiff(:, :, :),1)));
        title(['Same-Different Word Pairs Cosine Similarity for Block ', num2str(iBlock-1)])
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
        hold on
        plot(get(gca, 'xlim'), [timeZero timeZero], 'k', 'LineWidth', LT)
        plot([timeZero timeZero], get(gca, 'ylim'), 'k', 'LineWidth', LT)
        
        %%- Save Image
%         print(figureFile, '-dpng', '-r0')
%         savefig(figureFile)
        
        pause(0.1);
    end % loop through blocks
end % loop through sessions