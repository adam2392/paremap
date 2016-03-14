function responseTime = plotEventMetaData(events, subplotIndex, figMeta, NUMSUBPLOTS, TRIGGER)
    %% EXTRACT DATA FROM EVENTS
    %%%%%%%%%%%%%%%%%%%%%RECENTER THE TIMES WRT MSTIME %%%%%%%%%%%%%%%%
    eTimeOn = [events.mstime];
    fixationOn = ([events.fixationOnTime] - eTimeOn);     % fixation on time
    fixationOff = ([events.fixationOffTime] - eTimeOn);   % fixation off time
    matchOnTime = ([events.matchOnTime] - eTimeOn);
    probeOffTime = ([events.probeOffTime] - eTimeOn);

    % CONVERT ALL TIMES TO SECONDS
    eTimeOnS = (eTimeOn-eTimeOn)/1000;           % time on (seconds)
    fixationOnS = fixationOn/1000;     % fixation on time
    fixationOffS = fixationOff/1000;   % fixation off time
    matchOnTimeS = matchOnTime/1000;   
    probeOffTimeS = probeOffTime/1000;
    responseTimeS = [events.responseTime]/1000;

    %%- Create Histograms of events metadata
    labeloffset = 50;

    %%%%%% ONLY GET RESPONSE TIMES > 0 
    responseTimeS = responseTimeS(responseTimeS>0);
    
    % extract the information from a histogram function
    [n1, xout1] = hist(fixationOnS);
    [n2, xout2] = hist(eTimeOnS);
    [n3, xout3] = hist(fixationOffS);
    [n4, xout4] = hist(matchOnTimeS);
    [n5, xout5] = hist(probeOffTimeS);
    [n6, xout6] = hist(responseTimeS);
    
    %% BAR CHARTS
    figure(figMeta);
    %%-01 plot fixation on/off times - should be < 0
    subplot(NUMSUBPLOTS, 1,subplotIndex)
%     bar(xout1,n1,'r');          % fixation on
%     grid; hold on;
%     bar(xout3, n3, 'b');        % fixation off
%     bar(xout2,n2,'g'); grid;    % probeword on
    bar(xout6, n6, 'b'); grid;  % response times

%     set(gca, 'XLim', [min(xout1)-0.05, max(xout6)+0.05])
    set(gca, 'XLim', [0, max(xout6) + 0.05])
    
    % add text
%     text(mean(fixationOnS),max(n1)+labeloffset, 'fixation On')
%     text(mean(xout3), max(n3)+labeloffset, 'fixation off')
%     text(0, max(n2)+2*labeloffset,'probeWordOn')
    text(mean(xout6), max(n6)+labeloffset*6, 'response Time')
    title(sprintf('Plotting Events meta data for %s with %s events', TRIGGER, num2str(length(events))))


    responseTime = responseTimeS;
    % add legend
%     legend(...
%     'fixation On', ...
%     'fixation off', ...
%     'probewordOn', ...
%     'response Time', ...
%     'Location', 'northeast')
%     legend()
end