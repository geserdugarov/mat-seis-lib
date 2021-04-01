%% *mat2vec*
% Transform upper-right triangle of a square matrix into a row vector 

%%
% *Input:*

% Mat   - input square matrix
% shift - [scalar and positive] shift to the right from the main diagonal 

%%
% *Output:*

% Vec   - row vector

%%
% *Author:* Vladimir Grechka 2013

%%
function [Vec] = mat2vec(Mat, shift) 
%% Settings 
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
% narginTrue = nargin(thisFileName);   
% narginchk(narginTrue, narginTrue);

N = size(Mat,1);
if N ~= size(Mat,2)
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    display(Mat)
    fprintf('>>> Input matrix is not square \n \n');
      error('>>> STOP');
end;

if nargin < 2 
    shift = 0;
else
    if isempty(shift) == 1
        shift = 0;
    elseif shift < 0  ||  shift >= N
        fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
        fprintf(['>>> The second input variable ''shift'' (= %g) is supposed to obey', ...
                 ' the inequality 0 < shift <= %g \n \n'], shift, N);
          error('>>> STOP');
    end;
end;

%% Construct the output vector
i1 = 0;
for i=1 : N-shift
    i2 = i1 + length(i+shift : N);   
    Vec(1, i1+1 : i2) = Mat(i, i+shift : N);
    i1 = i2;
end;  

end    % of the function