% function varargout = get_filenames_spm5(varargin)
% 
% Function to get file names using SPM5. Supports both wildcard and
% regular expressions for file names (not for paths).
%
% INPUT:
%       file_path (1 x n string)
%       optional: filenames with wildcards (e.g. *.img or rf*.img), or 
%                 filenames with regular expressions starting with 
%                 REGEXP: (e.g. 'REGEXP:^rf.*\.img$')
%
% OUTPUT: filenames as n x m char array (n = number of files) with space padding at end
%
% Example calls:
%       fnames = get_filenames_spm5('/home/resultsdir','*.nii') 
%           will select all nifti files in the provided folder
%
%       fnames = get_filenames_spm5('/home/resultsdir/*.nii')
%           same result as above
%
%       fnames = get_filenames_spm5('/home/resultsdir','REGEXP:^rf.*\.img$')
%           will select all nifti files in the provided folder starting
%           with rf and ending with .img

function varargout = get_filenames_spm5(varargin)

if nargin == 1
    [fpath,fname,fext] = fileparts(varargin{1});
    varargin{1} = fpath;
    varargin{2} = [fname fext];
end

if length(varargin) >= 2 && strncmp(varargin{2},'REGEXP:',7)
    fname_regexp = varargin{2}(8:end);
else
    fname_regexp = wildcard2regexp(varargin{2});
end

varargout{1} = char(spm_select('Fplist',varargin{1},fname_regexp));