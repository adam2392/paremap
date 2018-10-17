function eventTrigger = timeLockEvents(events, timeLock)
    %%- Input to gete_ms
    %%- Dependent only on eventsTriggerXlim: These stay the same regardless of how we process events
    eventTrigger = events;

    % offset to synchronize with certain event trigger (e.g. vocalization of
    % word, or matchword coming on)
    if strcmp(timeLock, 'vocalization'),
        for iEvent=1:length(eventTrigger),
            eventTrigger(iEvent).mstime = eventTrigger(iEvent).mstime + eventTrigger(iEvent).responseTime;
            eventTrigger(iEvent).eegoffset = eventTrigger(iEvent).eegoffset + round(eventTrigger(iEvent).responseTime);
        end
        disp('Time Locking to VOCALIZATION');
    elseif strcmp(timeLock, 'matchword'),
        for iEvent=1:length(eventTrigger),
            eventTrigger(iEvent).eegoffset = eventTrigger(iEvent).eegoffset + round(eventTrigger(iEvent).matchOnTime - eventTrigger(iEvent).mstime);
            eventTrigger(iEvent).mstime = eventTrigger(iEvent).mstime + (eventTrigger(iEvent).matchOnTime - eventTrigger(iEvent).mstime);
        end
        disp('MATCHWORD');
    elseif strcmp(timeLock, 'probeword') % Settings for probewordon synchronization
        disp('PROBE WORD ON');
    else
        error('not set correctly.');
    end
end