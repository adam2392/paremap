%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script evcNetwork.m
%
% Description: Builds a graph based on the eigenvector centrality from the
% "4 assumed" regions in a seizure patient recording (ii, pre, sz, post).
% Computes a distance metric between each EVC in time to create 
%
% Input: the graph produce by 'evcNetwork.m' 
%
% Output: the different clusters that each EVC in time belongs to
%
% Author: Adam Li
%
% Ver.: 1.0 - Date: 08/21/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 01: Create Settings for Running Algorithm

addpath('/home/adamli/MATLAB/code_adam/burns_adapted/EVC/chr');
addpath('/home/adamli/MATLAB/code_adam/burns_adapted/EVC/pwr');
    
addpath('/Users/adam2392/Documents/MATLAB/Johns Hopkins/epilepsy/code_adam/burns_adapted/EVC/chr')
addpath('/Users/adam2392/Documents/MATLAB/Johns Hopkins/epilepsy/code_adam/burns_adapted/EVC/pwr')

Patient = 'PY04N007';
szindices = '61'; %[61 74 81];
mattype = 'q'; % q, or pwr_q
if strcmp(mattype, 'pwr_q')
    mat = 'pwr_'
else
    mat = ''
end
%% 02: Load in the EVC mat files 
%%%%%% Load in EVC files from /EVC/pwr or /EVC/chr
eval(['load ', mattype, szindices, 'pre_', Patient])    %%% Load in qpre
eval(['qpre = q',szindices, 'pre;'])  

eval(['load ', mattype, szindices, 'sz_', Patient])    %%% Load in qsz
eval(['qsz = q',szindices, 'sz;'])  

eval(['load ', mattype, szindices, 'post_', Patient])    %%% Load in qpost
eval(['qpost = q',szindices, 'post;'])  

eval(['load ', mattype, 'AllII_', Patient])    %%% Load in qAllII
try
    eval('qii = AllLeadSVD;')
catch error
    eval('qii = qAllII;')
end

% clear out unneeded variables after loading them in
clear(strcat(mattype, szindices, '*'), 'AllLeadSVD', 'qAllII') 
%% 03a: PreProcess EVC files and get certain time ranges 
%%% Timescale the EVC files to a certain length
rescaled = 1;

% first grab a smaller section of qii
qii = qii(300:799,:);

if rescaled
    offset = 150;
    
    % store duration of the region
    durpre = size(qpre, 1);        
    durpost = size(qpost, 1);
    dursz = size(qsz, 1);
    
    % create evenly spacedintervals 
    totallen = 500; % total length of the area of analysis (e.g. 500 seconds)
    intervalpre = linspace(1, durpre, offset);
    intervalpost = linspace(1, durpost, offset);
    intervalsz = linspace(1, dursz, totallen - 2*offset);   % create an interval of 500 points

    %%% Scale pre, seize, and post
    if dursz ~= totallen-offset %%%% Scale up and interpolate
       qsz = ((interp1(1:dursz, qsz, intervalsz, 'linear')));
    end
    if durpre ~= offset
       qpre = ((interp1(1:durpre, qpre, intervalpre, 'linear')));
    end
    if durpost ~= offset
       qpost = ((interp1(1:durpost, qpost, intervalpost, 'linear')));
    end    
end

%% 04: Create graph
% concatenate 
qcat = [qii; qpre; qsz; qpost];

%%%% graph will be a nxn matrix, where n is the length of qcat
D = pdist(qcat);
graph = squareform(D);

try
    clusters = spectralPartition(graph);
catch error
    
end