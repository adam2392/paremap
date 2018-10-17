function A = read_svd_array(filename,nch)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function: A = read_svd_array(filename,nch)
%
% Description: It accesses the file named "filename", reads data from it,
%              and format the data into a sequence of columns, each with
%              "nch" rows.
%             
% Input:    filename - Name of the file where the columns are stored. It
%                      must be a string of characters.
%
%           nch      - Number of rows of each column. It must be a positive
%                      integer.
%
% Output:   A - 2D array that stores the sequence of columns extracted from
%               the input file.
%			
%
%
% Author: S. Santaniello
%
% Ver.: 1.0 - Date: 11/16/2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%--------------------------------------------------------------------------
% check if the correct number and type of data has been passed
%--------------------------------------------------------------------------
if (nargin<2), error('Error: no enough input'); end

if (isempty(filename) || ~ischar(filename))
    error('Error: name of the source file not valid');
else
    % open the file
    fid = fopen(filename,'rb');
    if (fid<0)
        error('Error: file not opened correctly or not valid');
    end
end

if (~isreal(nch) || isempty(nch) || length(nch)>1 || nch<1)
    error('Error: number of rows not valid');
else
    nch = round(nch);
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% main loop
%--------------------------------------------------------------------------
fseek(fid,0,'bof');
A =fread(fid,[nch inf],'single');
fclose(fid);
%--------------------------------------------------------------------------
