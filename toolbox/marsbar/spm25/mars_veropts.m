function varargout = mars_veropts(arg, varargin)
% Returns SPM version specific parameters for MARSBAR with SPM25.
% SPM25 is the name for the development version of SPM12; it is presumed to function like SPM12.
% Therefore this is a minimal script which returns the appropriate configuration for the "SPM25" version.
% FORMAT varargout = mars_veropts(arg, varargin)

if nargin < 1
    varargout = {};
    return
end

switch lower(arg)
    case 'defaults'
        global defaults
        if isempty(defaults)
            defaults = spm('defaults','FMRI');
        end
        varargout = {defaults};
        
    case 'template_ext' 
        varargout = {'.nii'}; % NIfTI format
        
    case 'pref_img_out_ext'
        varargout = {'.nii'}; % Preferred output format
        
    case 'des_conf'
        varargout = {'SPM.mat'}; % SPM design file
        
    case 'design_filter_spec'
        varargout = {{'SPM.mat','SPM.mat'}}; % Design filter spec
        
    case 'flip_option'
        varargout = {0}; % No flip needed for NIfTI
        
    otherwise
        error(['Unsupported option: ' arg]);
end