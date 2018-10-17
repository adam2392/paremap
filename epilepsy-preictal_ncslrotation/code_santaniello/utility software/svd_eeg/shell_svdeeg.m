function shell_svdeeg(pathval,filename,pos)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function: shell_svdeeg(pathval,filename,pos)
%
% Description: This function is a shell-interface for data processing.
%              Input arguments (char []*) are used to locate the source
%              file of multi-channel iEEG recordings, and are passed to the
%              routines that compute the connectivity matrices and the
%              correspondent singular value decomposition. Results are
%              stored in *.dat files.
%             
% Input:    pathval      - Path to the directory where the destination and
%                          source files are stored. It must be a string of
%                          characters.
%
%           filename     - Name of the *.rec file where the iEEG recordings
%                          are extracted from. It must be a string of
%                          characters.
%
%           pos          - pointer to the last second of data accessed
%                          (counted from 0), i.e., the first byte to-be-
%                          extracted is in position: pos x num_channels x
%                          #-samples-per-second x #-bytes-per-sample + 1.
%                          It is an ASCII string that must be converted to
%                          integer.
%
% Output:   no output returned.
%			
%
%
% Author: S. Santaniello
%
% Ver.: 5.0 - Date: 11/05/2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%--------------------------------------------------------------------------
% check if the correct number and type of data has been passed
%--------------------------------------------------------------------------
if (nargin<3), error('Error: no enough input'); end

if (isempty(pathval) || ~ischar(pathval) || ~isdir(pathval))
    error('Error: path is not a valid directory');
end

if (isempty(filename) || ~ischar(filename) || ~exist(sprintf('%s/%s.rec',pathval,filename),'file'))
    error('Error: name of the source file not valid or file not found');
end

if (isempty(pos) || ~ischar(pos))
    error('Error: position not valid');
else
    pos = round(str2double(pos));
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% shell
%--------------------------------------------------------------------------
% extract the formatting information about the iEEG recordings
hdr = readhdr(sprintf('%s/eeg.hdr',pathval));

% extract the number of recording channels
num_channels = hdr.file_fmt.Numb_chans; 

% extract the format of the data in the *.rec file
format_file  = hdr.file_fmt.File_format;    

% extract the offset in the *.rec file required for the EDF header
offset  = hdr.file_fmt.Data_offset;    

% extract the sampling frequency (in # of samples) in the *.rec file
Fs = hdr.file_fmt.Samp_rate;

% compute the cross-power and coherence based connectivity matrices
power_coherence(pathval,filename,format_file,num_channels,Fs,offset,pos);

% set a pointer to the channels that must be used for singular value
% decomposition
subj_name = pathval(strfind(pathval,'PY'):strfind(pathval,'PY')+7);
switch subj_name
    case 'PY04N007'
        included_chn = [1:79 81:87];
    
    case 'PY04N008'
        included_chn = 1:40;
    
    case 'PY04N009'
        if (~isempty(strfind(pathval,'edf1')))
            included_chn = [1:2 5:24 27:59 63:73];
        else
            included_chn = [1:2 5:24 27:59 62:72];
        end
        
    case 'PY04N012'
        if (~isempty(strfind(pathval,'edf0')))
            included_chn = [1:34 51:98];
        else if (~isempty(strfind(pathval,'edf2')))
                included_chn = 1:82;
            else
                included_chn = 1:num_channels;
            end
        end
        
    case 'PY04N013'
        included_chn = [1:10 12:61 63 65 69:89];
    
    case 'PY04N015'
        included_chn = [1:2 5:83 85:95];
    
    case 'PY05N004'
        if (~isempty(strfind(pathval,'edf0')))
            included_chn = [1:12 14 17:37 39:63 65:82 84:100];
        else
            included_chn = [1:12 14 17:37 39:81 83:99];
        end
        
    case 'PY05N005'
        included_chn = [1:21 23:38 41:83 89:118];
    
    case 'PY05N006'
        included_chn = [1:2 5:24];
    
    case 'PY05N007'
        included_chn = [1:2 5:7 10:14 17:48];
    
    case 'PY05N019'
        included_chn = [1:97 100:115];
    
    case 'PY11N003'
        included_chn = 1:117;
    
    case 'PY11N004'
        included_chn = [2 4:15 17:31 33:106 108:117];
    
    case 'PY11N005'
        if (~isempty(strfind(pathval,'edf0')))
            included_chn = [1:11 13:33 35:74];
        else
            included_chn = [1:32 34:73];
        end
    
    case 'PY11N006'
        included_chn = [1:16 18:66];
    
    case 'PY11N007'
        included_chn = [1:5 7:18 20:33 35:48];
    
    case 'PY11N008'
        included_chn = [1:113 115:126];
    
    case 'PY11N009'
        included_chn = [1:22 24:30 32:81 85:90 92:105];
    
    case 'PY11N010'
        included_chn = [1:77 79:85 87:88];
    
    case 'PY11N011'
        included_chn = 6:32;
    
    case 'PY11N012'
        included_chn = [1:55 57:89 91:110];
    
    case 'PY11N013'
        included_chn = [1:111 113:119];
    
    case 'PY11N014'
        included_chn = [1:21 23:40 42:46 48:52 55:68 70:81 83:111];
    
    case 'PY11N015'
        included_chn = [5:15 17:32];
end

% singular value decomposition
svd_decomposition(pathval,filename,num_channels,included_chn); 
%--------------------------------------------------------------------------
