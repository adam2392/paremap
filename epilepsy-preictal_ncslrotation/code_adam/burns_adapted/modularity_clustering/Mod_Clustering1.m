%filename is the name of the patient file in quotations
%seiz is a 1xn row vector containing the numbers of n seizures to be
%processed (i.e. [1 3 4] would process the first, third, and fourth
%seizures.

%It was easier to use cells and index them for each seizure (i.e. ranks{1}
%is the first seizure you selected, ranks{2} is the second, etc.). This
%way, the number of seizures can be varied and it's easier to create loops

function [C,all_events,clusters,correlationsB,avg_rank_signals,mod] = Mod_Clustering1(patientnum,seiz,measure,freqband)

%Load Files
if strcmp(freqband,'beta')
    band = '13-25';
elseif strcmp(freqband,'gamma')
    band = '62-100';
end
    
load(sprintf(['pat_%d_svd_seizure_blocks' band],patientnum));
load(sprintf('pat_%d_g_electrodes',patientnum))

r = size(allFirstSingVecs,1);
w = [];

%Create row vector w with each entry q the length of a chunk
q = 1;
for i = seiz
    w(q) = (end_marks(i)+120) - (start_marks(i)-120);
    q = q+1;
end

%Preallocate the size of the ranks matrices
ranks = [];
for i = 1:size(w,2)
    ranks{i} = ones(w(i),r);
end

%Inserts actual entries in ranks matrices
B = [];
rank = [];
pos = 1;
for k = seiz %loop through each seizure
    if end_marks(k)+120 > block_marks(k+1)
        P = block_marks(k+1);
    else
        P = end_marks(k)+120;
    end

    for i = (start_marks(k)-120):P
        [B{pos},rank{pos}] = sort(allFirstSingVecs(:,i),'descend');
        for j = 1:r
            ranks{pos}((i-(start_marks(k)-120)+1),j) = find(rank{pos} == j);
        end
    end
    pos = pos + 1;
end

%DC OFFSET OF EACH
ranks_smooth = [];
for i = 1:size(ranks,2)
    for j = 1:r;
        ranks_smooth{i}(:,j) = (ranks{i}(:,j)) - mean(ranks{i}(:,j));
    end
end



%CONCATENATION
C = [];
for i = 1:size(ranks_smooth,2)
    C = cat(1,C,ranks_smooth{i});
end

%SMOOTHING OUT
for i = 1:r;
    C(:,i) = smoothout2(C(:,i)); %- mean(smoothout2(C(:,i)));
end

%MAKING ALL_EVENTS; could have been in the same loop as concatenation but
%I made it separate for clarity
all_events = ranks{1};
for i = 2:size(ranks,2)
    all_events = cat(1,all_events,ranks{i});
end

%C is the matrix of rank signals- with the channels as COLUMNS and time in
%seconds as ROWS (you may have to take the transpose of this matrix,
%depending on what you are doing with the data)-
%after the average has been subtracted and it has been smoothed out

%all_events is the raw ranked signal of all events (before the average was
%subtracted and it was smoothed and stuff) (also channels as COLUMNS and
%time as ROWS

%Other files: matplot, run_KL, the other KL thing, smoothout2, Reorganize

%MAKING THE CORRELATION MATRIX; for type input, either 'mse' or 'xcorr'
[w1,r] = size(all_events);
correlationsB = ones(r,r);

if strcmp(measure,'mse') == 1
    for i = 1:r
        for j = 1:r
            correlationsB(i,j) = -mean((C(:,i) - C(:,j)).^2);
        end
    end
else
    for i = 1:r
        for j = 1:r
            correlationsB(i,j) = max(xcorr(C(:,i),C(:,j)));
        end
    end
end

i = 1;
j = 1;
for k = 1:r
    correlationsB(i,j) = 1;
    i = i + 1;
    j = j + 1;
end

% % % %PLOTTING UNSHUFFLED CORRELATION MATRIX- if you want to see what it looks
% % % %like without being shuffled- you can comment this out if you want
matplot(correlationsB)
pause

%Run KL to get to 2^n clusters and plot them.
%Austin 7/17/2013 - I changed this code from being 16 all the time to being
%able to vary the number of clusters being created based on how many
%electrodes there are, since the new code allows for the removal of
%electrodes and the new number might not allow for 16 clusters

COMTY = cluster_jl(correlationsB);
commID = COMTY.COM{1};
comms = {};
comms{size(COMTY.SIZE{1},2)} = [];

for i = 1:size(commID,2)
    comms{commID(i)} = [comms{commID(i)} i];
end

maximum = 0;
for i = 1:size(comms,2)
    if size(comms{i},2) > maximum
        maximum = size(comms{i},2);
    end
end

clusters = zeros(size(comms,2),maximum);
for i = 1:size(comms,2)
    clusters(i,1:size(comms{i},2)) = comms{i};
end

mod = COMTY.MOD(1);

%Smoothout all_events before plotting it
for i = 1:size(all_events,2)
    all_events(:,i) = smoothout3(all_events(:,i));
end

%Pull all the rows relevant to each cluster, average them, and store all
%the averages in a new matrix
avg_rank_signals = []; %- the new matrix to store it in
[len, wid] = size(clusters);
for i = 1:len
    indices = clusters(i,:); %- each row of clusters
    indices(indices == 0) = []; %- not counting any of the 0's in any given row
    cluster = [];
    for j = 1:length(indices)
        cluster(j,:) = (all_events(:,indices(j))); %- Pull the relevant columns of all_events for the electrodes in each cluster
    end
    avg_rank_signals(i,:) = smoothout3(mean(cluster,1)); %- average and smoothout each cluster into one average rank signal and store that in the new matrix
end

%Create an array containing the lengths of each seizure
q = 1;
seizlen = zeros(1,size(seiz,2));
for i = seiz
    seizlen(q) = end_marks(i)-start_marks(i);
    q = q + 1;
end

%Preallocate the size of the arrays that will contain the start and end
%points of the seizures
starts = zeros(1,size(seiz,2));
ends = zeros(1,size(seiz,2));

%Enter start and ends points of seizures into separate arrays
q = 1;
for i = seiz
    mostseiz = 0;
    allseiz = 0;
    for j = 1:q
        if j == q
            allseiz = allseiz + seizlen(j);
        else
            allseiz = allseiz + seizlen(j);
            mostseiz = mostseiz + seizlen(j);
        end
    end            
    starts(q) = 120*(2.*q-1) + mostseiz;
    ends(q) = 120*(2.*q-1) + allseiz;
    q = q + 1;
end

%Plot the average rank signatures of all clusters
t = 1:1:w1;
figure(2);
plot(t,avg_rank_signals)
hold on;
title('Average Rank Signal of All Clusters')

%Plotting lines to mark the beginning and end of each seizure; green lines
%mark the beginning, red lines mark the end, black lines are the borders
%between recordings
for i = 1:size(seiz,2)
    line([starts(i) starts(i)],[ylim],'color','g','linewidth',1.5)
    line([ends(i) ends(i)],[ylim],'color','r','linewidth',1.5)
    if i == size(seiz,2)
        continue
    else
        line([(ends(i)+120) (ends(i)+120)],[ylim],'color','k','linewidth',1.5)
    end
end
hold off

% g_electrodes
clusters
starts
ends
all_events

end
% Summary of what you get out of this:

% Clusters: a matrix where each row is a cluster. The row indices from the
% original channels that go in each cluster are in each row. Zeros take up
% the extra space in each row because not all clusters are the same size

% correlationsB: the "correlation matrix" that has the peak correlation
% between each pair of signals- unshuffled