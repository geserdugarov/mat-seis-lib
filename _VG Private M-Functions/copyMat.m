%% *copyMat*
% Copy elements of one matrix onto another

%%
% *Input:*

% mat1, mat2 - [:, :] two input matrixes

%%
% *Output:*

% mat3       - [:, :] matrix containing the elements of mat1 copied onto mat2; 
%              size(mat3) = size(mat2)

%%
% *Author:* Vladimir Grechka 2012 - 2014

%%
function [mat3] = copyMat(mat1, mat2)
%% Settings  
[~, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

%% Map mat2 onto mat3
mat3 = mat2;

%% Map mat1 onto mat3
n1 = size(mat1);   n2 = size(mat2);
min1 = min([n1(1), n2(1)]);  
min2 = min([n1(2), n2(2)]);
mat3(1:min1, 1:min2) = mat1(1:min1, 1:min2);

% %% Loop over the smallest dimensions of the input matrices
% for i1 = 1:min([n1(1), n2(1)])
%     for i2 = 1:min([n1(2), n2(2)])
%         mat3(i1,i2) = mat1(i1,i2);
%     end;
% end;

end   % of the function