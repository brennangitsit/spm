% function output = transres_confusion_matrix(decoding_out, chancelevel, varargin)
% 
% Get a matrix of how often each label has been confused by other labels.
% This function will calculate the per-class accuracy, i.e. balanced for
% the number of occurrences of each label. The output will be percentage 
% of appearance of the class, i.e. each row will sum up to 100. To change
% this, modify the commented line in the code and save under a different
% name.
%
% The output will be an NxN matrix where n is the number of unique labels.
% The columns will represent the predicted labels, whereas the rows will
% represent the true labels. The output is sorted by label number,
% from low to high.
%
% To use this transformation, use 
%
%   cfg.results.output = {'confusion_matrix'}
%
% NOTE: Many classifiers, including libsvm, will in multiclass settings
% assign samples for which no unique label can be predicted to one class 
% (libsvm eg to the first class). This leads to a potentially heavy bias 
% of the first class. To take this into account, use
%
%   cfg.results.output = {'confusion_matrix_plus_undecided'}
%
% which will create one additional 'undecided' column to the confusion 
% matrix that contains the number of all samples for which no unique
% prediction was possible.
%
% Martin, 2014-04-23
%
% See also decoding_transform_results transres_accuracy_matrix transres_accuracy_pairwise

% Update KG 2021-03-16: added confusion_matrix_plus_undecided to header
% Update MH 2017-03-15:
% Error in description: flipped columns and rows, now corrected

function output = transres_confusion_matrix(decoding_out, chancelevel, varargin)

predicted_labels =  vertcat(decoding_out.predicted_labels);
true_labels = vertcat(decoding_out.true_labels);

labels = uniqueq(true_labels);
n_labels = size(labels,1);

output = zeros(n_labels); % init

for i_label = 1:n_labels
    labelfilt = true_labels == labels(i_label);
    curr_predicted_labels = predicted_labels(labelfilt);
    curr_n_true = sum(labelfilt);
    for j_label = 1:n_labels
        output(i_label,j_label) = 100 * (1/curr_n_true) * sum(curr_predicted_labels==labels(j_label)); % deactivate me if you want raw results
        % output(i_label,j_label) = sum(curr_predicted_labels==labels(j_label)); % activate me if you want raw results
    end
end
output = {output};