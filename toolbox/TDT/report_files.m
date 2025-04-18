% function report_files(cfg,n_steps,inputfilenames_fid)
% 
% This function prints file names for a certain decoding to the screen and
% writes them to a given file if requested. It is an integral part of the
% decoding toolbox and should not be called directly.

function report_files(cfg,n_steps,inputfilenames_fid)
global verbose
local_verbose = 2; % verbosity of current file

if ~cfg.results.write
    if verbose < local_verbose
        dispv(1, ['report_files:Skipping reporting files because cfg.results.write=1 and (global verbose) ' num2str(verbose) ' >= ' num2str(local_verbose) ' (local verbosity level)'])
        return % nothing to do
    end
    inputfilenames_fid = ''; % this will prevent writing to file/only to screen
end

% Find common string in all files to print this only once
fnames = cfg.files.name;
fnames_char = char(fnames);
n_str = size(fnames_char,2); % maximum string length
for i_str = 1:n_str
    match = strncmp(fnames{1},fnames(2:end),i_str);
    if ~all(match)
        n_match = i_str-1;
        break
    end
end
filestart = fnames_char(1,1:n_match); % common file start
filerest =  fnames_char(:, n_match+1:end); % get not common part

for i_step = 1:n_steps

    % Get indices for training
    i_train = find(cfg.design.train(:, i_step) > 0);
    % Get indices for testing
    i_test = find(cfg.design.test(:, i_step) > 0);

    if isfield(cfg, 'sn')
        text = sprintf('Subject %i, Decoding Nr %i', cfg.sn, i_step);
    else
        text = sprintf('Decoding Nr %i', i_step);
    end
    dispv(local_verbose, '%s', text)
    fprintf(inputfilenames_fid, '%s\n', text);

    if n_match > 0
        cont = '...';
        text = sprintf('  File Start: %s%s\n', filestart, cont);
        dispv(local_verbose, '%s', text)
        fprintf(inputfilenames_fid, '%s\n', text);
    else
        cont = '';
    end

    for curr_i_train = i_train'
        if iscell(cfg.design.label(curr_i_train, i_step)) % for multitarget
            text = sprintf('  File Train (%s): %s%s', num2str(cfg.design.label{curr_i_train, i_step}), cont, filerest(curr_i_train,:));    
        else
            text = sprintf('  File Train %i: %s%s', cfg.design.label(curr_i_train, i_step), cont, filerest(curr_i_train,:));
        end
        dispv(local_verbose, '%s', text)
        fprintf(inputfilenames_fid, '%s\n', text);
    end
    fprintf(inputfilenames_fid, '\n');


    for curr_i_test = i_test'
        if iscell(cfg.design.label(curr_i_test, i_step)) % for multitarget
            text = sprintf('  File Test (%s): %s%s', num2str(cfg.design.label{curr_i_test, i_step}), cont, filerest(curr_i_test,:));   
        else
            text = sprintf('  File Test %i: %s%s', cfg.design.label(curr_i_test, i_step), cont, filerest(curr_i_test,:));
        end
        dispv(local_verbose, '%s', text)
        fprintf(inputfilenames_fid, '%s\n', text);
    end
    fprintf(inputfilenames_fid, '\n');

end