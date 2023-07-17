function [ out ] = computeXTX( X, nblocks )
% NEWFUN
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
% Optional
%--------------------------------------------------------------------------
% OUTPUT
% 
%--------------------------------------------------------------------------
% EXAMPLES
% 
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------

%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist( 'opt1', 'var' )
   % Default value
   opt1 = 0;
end

%%  Main Function Loop
%--------------------------------------------------------------------------
n = size(X,1);
nsubj_per_block = floor(n/nblocks);
nsubsinlaskblock = n - (nsubj_per_block)*(nblocks-1);

XTX = zeros(n,n);

% Initialize the indices for the blocks
block_subject_indices = (1:nsubj_per_block) - nsubj_per_block;

for I = 1:nblocks-1
    block_subject_indices = block_subject_indices + nsubj_per_block;
    XTX = XTX + X(block_subject_indices, :)'*X(block_subject_indices, :);
end



end

