%% *removeNaNs*
% Remove rows or columns containing NaNs from a matrix 

%%
% *Input:*

% matIn  - input matrix
%          (*) ndims(matIn) is assumed to be equal to 2 
% dim    - dimension along which NaNs are to be removed:
%          dim = 1 for columns and dim = 2 for rows

%%
% *Output:*

% matOut - output matrix from which the rows or columns containing NaNs are removed

%%
% *Author:* Vladimir Grechka 2013 - 2014

%%
function [matOut] = removeNaNs(matIn, dim)
%% Settings
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

if ismatrix(matIn) == 0
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Input ''matIn'' is not a matrix: ndims(matIn) = %g \n', ndims(matIn));
      error('>>> STOP');
end;

matOut = matIn;

%% Remove NaNs
if dim == 1
    matOut(:, any(isnan(matOut),1)) = [];
elseif dim == 2
    matOut(any(isnan(matOut),2), :) = [];
else    
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Input parameter dim = %g. Its legitimate values are dim = 1 or 2 \n', dim);
      error('>>> STOP');
end;
matOut = squeeze(matOut);

end  % of the function