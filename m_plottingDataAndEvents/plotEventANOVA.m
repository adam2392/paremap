function plotEventANOVA(anovaPowMat, binnedWaveT, freqBandAr, titleStr, ...
    figAnova)
    %% EXTRACT DATA FROM EVENTS
    freqBandYtickLabels = {freqBandAr.label};
    %% Actually Run ANOVA For This Channel
    spectMat = anovaPowMat{1};
    
    anovaMat = zeros(size(spectMat,2), size(spectMat,3));
    for freq=1:size(spectMat,2)
        for time=1:size(spectMat,3)
            y = [];
            groups = [];
            for i=1:length(anovaPowMat)
                % create vector of events we want to test
                y = [y; anovaPowMat{i}(:,freq,time)];
                
                % set groups
                group = ones(size(anovaPowMat{i}(:,freq,time)))*i;
                groups = [groups; group];
            end
            
            % compute p-value for ANOVA
            p = anovan(y, groups, 'display','off');
            %%%% 2nd group with target words,...
            
            % add to ANOVA matrix of p-values
            anovaMat(freq,time) = p;
        end
    end
    
    anovaMat(anovaMat > 0.05) = 0.5;
    min(anovaMat(:))
    [r,c] = find(anovaMat == min(anovaMat(:)),1)
    %% Spect Map of ANOVA
    figure(figAnova);
%     subplot(NUMSUBPLOTS, 1,subplotIndex)  
    hImg    = imagesc(binnedWaveT,[1:7],anovaMat); 
    hold on;  colormap(jet);

    hCbar = colorbar('east');
    set(hCbar,'ycolor',[1 1 1]*0.1, 'YAxisLocation', 'right')
            
    % set the heat map settings
    set(gca,'ytick',[1:7],'yticklabel',freqBandYtickLabels)
    set(gca,'tickdir','out','YDir','normal'); % spectrogram should have low freq on the bottom
    title(sprintf('Plotting ANOVA: %s', titleStr))

    % label the ranges for each band
    xlim=get(gca, 'xlim');
    rangeF = [1:7; 2:8]';
    for iFB=1:length(freqBandAr),
        hB = plot(xlim, [1 1]*rangeF(iFB,1)+0.5, '--', 'linewidth', 2, 'color', [1 1 1]*.7);
    end
end