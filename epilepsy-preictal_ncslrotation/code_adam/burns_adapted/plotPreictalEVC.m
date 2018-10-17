%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script plotPreictalEVC.m
%
% Description: plots the EVC's of preictal state
%
% Output:   no output returned.
%
% Author: Adam Li
%
% Ver.: 1.0 - Date: 09/22/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear
%% 01: Create Settings 
Patient = 'PY04N013';
FreqBands = {'alpha', 'beta', 'gamma', 'high'};
FreqBand = 'beta';
timeWindow = '7200'; 
% szindices = {'61' '74' '81'};
% szindices = {'124A', '124B', '125'}; %PY05N004
% szindices = {'762' '767' '777' '780' '783' '801'}; %PY04N013
% szindices = {'8', '20'}; %PY04N008 
concatenate = 1;
szindices = {'767' '777' '780' '783' '801'};

distanceFunc = 'cosine'; %or cityblock, correlation, sqeuclidean
options = statset('MaxIter', 3000);

addpath('/home/adamli/MATLAB/code_adam/burns_adapted/focus_electrode_analysis/7200');
addpath('/Users/adam2392/Documents/MATLAB/Johns Hopkins/epilepsy/code_adam/burns_adapted/focus_electrode_analysis/7200');

% loop through each frequency band
for f=1:length(FreqBands)
%     FreqBand = FreqBands{f}
    FreqBand = 'high';
    
    addpath(strcat('/home/adamli/MATLAB/code_adam/burns_adapted/focus_electrode_analysis/', FreqBand))
    addpath(strcat('/Users/adam2392/Documents/MATLAB/Johns Hopkins/epilepsy/code_adam/burns_adapted/focus_electrode_analysis/', FreqBand));

    % loop through each seizure event
    for i=2:length(szindices)
        seizurenum = szindices{i};

        % Load in the EVC data from 7200 out
        disp(['q' seizurenum 'pre_' Patient FreqBand])
        eval(['load ' 'q' seizurenum 'pre_' Patient FreqBand ';'])
        eval(['data = q' seizurenum 'pre_' Patient FreqBand ';'])  
        
        % load in EVC data of 600 seconds out
        eval(['load q' seizurenum 'pre_' Patient '_600;'])
        
        %% Grab a period of 600 seconds from each time Window
        % 7200 seconds out
        qpre7200 = data(1:600,:);
        
        % 3600 seconds out
        qpre3600 = data(4001:4600,:);
        
        % 1800 seconds out
        qpre1800 = data(end-1800+1:end-1800+600,:);
        
        % original 600 seconds out
        eval(['qpre600 = q' seizurenum 'pre;'])
        
        %% SVD and/or cluster using unsupervised on concatenated EVC
        if concatenate
            qallpre = [qpre7200;qpre3600;qpre1800;qpre600];

            [idxall, centsall] = kmeans(qallpre, 5, 'Distance', distanceFunc, 'Options', options);
            distancesall = pdist(centsall, 'cosine');
            distancesall = squareform(distancesall);

            save(strcat(Patient, '_', seizurenum, '_concatenateEVC_allpredistances'), 'distancesall', 'centsall', 'idxall', 'qallpre')

            Z = figure(2*length(szindices)+i);
            uniqueidx = unique(idxall);
            % Show firing rate for each preictal state
            for k=1:length(uniqueidx)
               showstate = uniqueidx(k);
               spiketrain = idxall == showstate;

               % compute the firing rate
               step = 5;
               for l=1:(length(spiketrain)-1)/step
                   if l==1
                       index = 1;
                   else
                       index = l*step;
                   end
                   showfiring(l) = sum(spiketrain(index:index+step))/step;
               end

               z = subplot(3, 2, k);
               plot(showfiring);
               title(['Firing rate for state: ' num2str(uniqueidx(k)) ' and seizure: ' seizurenum])
            end

            saveas(Z, ['allIDXFiringRate_' seizurenum '.bmp'])

            X = figure(21*i);
            plot(idxall);
            title(['Clusters for the ' num2str(i) 'th seizure all time windows concatenated'])
            xlabel('time (seconds)')
            ylabel('states')

            %% plotting distance matrices
            Y = figure(17*i);
            imagesc(distancesall);
            colorbar;
            title(['Cosine Distance Matrices for the ' num2str(i) 'th seizure all time windows concatenated'])
            
            % save dist matrix, cluster plot and firing rate
            currentdir = pwd;
            saveas(Y, [Patient '_' num2str(i) 'seizure_distancemat.jpg']);
            saveas(X, [Patient '_' num2str(i) 'seizure_clusterplot.jpg']);
            saveas(Z, [Patient '_' num2str(i) 'seizure_firingrates.jpg']);
            
            %% SVD
%         [dU,dS,dV] = svd(qallpre);
%         eigenvalues = diag(dS);
%         
%         E = figure(13*i);
%         plot(eigenvalues);
%         title(['Plot of eigenvalues for ' num2str(i) 'th seizure'])
%         
%         % pca analysis
%         coeff = pca(qallpre);
%         %%% Build a KD-tree to perform nearest neighbor search based on the
%         %%% Unsupervised PreCentroids
%         mdl = createns(coeff, 'NSMethod', 'kdtree', 'Distance', 'minkowski', 'BucketSize', 80);
%         clustersfrompca = knnsearch(mdl, qallpre, 'K', 1);
%         
%         F = figure(14*i);
%         plot(clustersfrompca)
%         title(['Plot of clusters using PCA'])
        
        %% generate cluster on each window separately using unsupervised k-means
        else
            qclusters = 5;
            %7200
            [idx7200, cents7200] = kmeans(qpre7200, qclusters, 'Distance', distanceFunc, 'Options', options);  
            distances7200 = pdist(cents7200, 'cosine');
            distances7200 = squareform(distances7200);

            %3600
            [idx3600, cents3600] = kmeans(qpre3600, qclusters, 'Distance', distanceFunc, 'Options', options);  
            distances3600 = pdist(cents3600, 'cosine');
            distances3600 = squareform(distances3600);

            %1800
            [idx1800, cents1800] = kmeans(qpre1800, qclusters, 'Distance', distanceFunc, 'Options', options);  
            distances1800 = pdist(cents1800, 'cosine');
            distances1800 = squareform(distances1800);

            %600
            [idx600, cents600] = kmeans(qpre600, qclusters, 'Distance', distanceFunc, 'Options', options);  
            distances600 = pdist(cents600, 'cosine');
            distances600 = squareform(distances600);

            allcents = [cents7200;cents3600;cents1800;cents600];
            overalldistance = pdist(allcents, 'cosine');
            overalldistance = squareform(overalldistance);
            save(strcat(Patient, '_', seizurenum, '_overallpredistances'), 'overalldistance', 'allcents')

    %         %%% Build a KD-tree to perform nearest neighbor search based on the
    %         %%% Unsupervised PreCentroids
    %         mdl = createns(precents, 'NSMethod', 'kdtree', 'Distance', 'minkowski', 'BucketSize', 80);
    %         newprecluster = knnsearch(mdl, predata, 'K', 1);
    %         unique(newprecluster);

            %% Plotting the EVC Heat Map
            A = figure(i);
            a = subplot(4, 1, 1);
            imagesc(transpose(qpre7200));
            title(['For the ' num2str(i) 'th seizure 7200 seconds before seizure'])
            colorbar;

            b = subplot(4, 1, 2);
            imagesc(transpose(qpre3600));
            title(['For the ' num2str(i) 'th seizure 3600 seconds before seizure'])
            colorbar;

            c = subplot(4, 1, 3);
            imagesc(transpose(qpre1800));
            title(['For the ' num2str(i) 'th seizure 1800 seconds before seizure'])
            colorbar;

            d = subplot(4, 1, 4);
            imagesc(transpose(qpre600));
            title(['For the ' num2str(i) 'th seizure 600 seconds before seizure'])
            colorbar;

            %% Plotting k-mean clusters
            B = figure(10*i);
            a = subplot(4, 1, 1);
            plot(idx7200);
            title(['Clusters for the ' num2str(i) 'th seizure 7200 seconds before seizure'])

            b = subplot(4, 1, 2);
            plot(idx3600);
            title(['Clusters for the ' num2str(i) 'th seizure 3600 seconds before seizure'])

            c = subplot(4, 1, 3);
            plot(idx1800);
            title(['Clusters for the ' num2str(i) 'th seizure 1800 seconds before seizure'])

            d = subplot(4, 1, 4);
            plot(idx600);
            title(['Clusters for the ' num2str(i) 'th seizure 600 seconds before seizure'])

            %% plotting distance matrices
    %         C = figure(11*i);
    %         title(['Distance Matrices for the ' num2str(i) 'th seizure'])
    %         a = subplot(2, 2, 1);
    %         imagesc(distances7200);
    %         colorbar;
    %         title('distance matrix 7200 seconds from seizure')
    %         
    %         b = subplot(2, 2, 2);
    %         imagesc(distances3600);
    %         colorbar;
    %         title('distance matrix 3600 seconds from seizure')
    %         
    %         c = subplot(2, 2, 3);
    %         imagesc(distances1800);
    %         colorbar;
    %         title('distance matrix 1800 seconds from seizure')
    %         
    %         d = subplot(2, 2, 4);
    %         imagesc(distances600);
    %         colorbar;
    %         title('distance matrix 600 seconds from seizure')

            % overall distance matrix
            D = figure(12*i);
            imagesc(overalldistance)
            title('Overall distance matrix')
            colorbar;
        end
    end
end


