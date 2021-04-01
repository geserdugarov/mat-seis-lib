%% *vec2mat*
% Transform a vector, which supposedly represents the upper-right triangle of a square symmetric 
% matrix, into the matrix 

%%
% *Input:*

% Vec - input vector
%       (*) length(Vec) should be equal to N*(N+1)/2, where integer N is the size of output 
%           symmetric matrix 

%%
% *Output:*

% Mat - symmetric matrix

%%
% *Author:* Vladimir Grechka 2012 

%%
function [Mat] = vec2mat(Vec)
%% Settings 
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

% Solve the equation length(Vec) = N*(N+1)/2 for N to verify that vector Vec has a legitimate 
% length 
D = 1 + 8*length(Vec);
N = (sqrt(D) - 1)/2; 
if isequal(floor(N), ceil(N)) == 0
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> The length %g of input vector produces a non-existent square matrix \n', ...
            length(Vec));
    fprintf('    that has non-integer size [%g, %g] \n \n', [N, N]);
      error('>>> STOP');
end;
 
%% Build the output matrix
Mat = zeros(N, N);
i1 = 0;
for i = 1:N
    i2 = i1 + length(i:N);   
    Mat(i, i:N) = Vec(i1+1:i2);
    i1 = i2;
end;  
Mat = symMat(Mat);

end    % of the function