%% *enhanceDiagDom*
% Sort columns of a matrix to enhance its diagonal dominance

%%
% *Input:*

% Mat   - input matrix 

%%
% *Output:*

% index - [1, size(Mat,2)] array of the indexes of sorted columns

%%
% *Author:* Vladimir Grechka 2012 

%%
function [index] = enhanceDiagDom(Mat)
%% Settings  
[~, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

%% Find the greatest absMat in each column in the order of absMat 
absMat = abs(Mat);
index  = zeros(1, size(Mat,2));

for i = 1:size(Mat,2)
    [row, col] = find(absMat == max(max(absMat)));
    index(1,col(1)) = row(1);

    % Remove the row [row(1),:] and column [:,col(1)] from absMat
    absMat(row(1),:) = NaN;   
    absMat(:,col(1)) = NaN;
end;

end    % of the function
