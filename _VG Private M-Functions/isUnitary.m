%% *isUnitary*
% Check whether input matrix is unitary within given tolerance

%%
% *Input:*

% Mat     - matrix
% tol     - tolerance of the proximity of R to a unitary matrix
% verbose - [scalar] equal to 1 or 0 to allow (verbose = 1) or suppress (verbose = 0) print out

%%
% *Output:*

% flag    - [scalar] equal to 1 if Mat is unitary within given tolerance and to 0 otherwise
% diffU   - [size(Mat)] matrix equal to Mat'*Mat - eye(size(Mat)) and quantifying the deviation 
%           of Mat from a unitary matrix

%%
% *Author:* Vladimir Grechka 2012 2014

%%
function [flag, diffU] = isUnitary(Mat, tol, verbose)  
%% Settings 
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
% narginTrue = nargin(thisFileName);   
% narginchk(narginTrue, narginTrue);

if nargin < 2;   tol = 1.e-6;   end;   % default for 'tol'
if nargin < 3;   verbose = 1;   end;   % default for 'verbose'

% Check whether Mat is square
[n, m] = size(Mat);
if n ~= m
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Input matrix is not square \n \n')
      error('>>> STOP');
end; 

%% Verify the unitarity
flag = 1;
diffU = Mat'*Mat - eye(size(Mat));
if max(max(abs(diffU))) > tol 
    flag = 0;
    if verbose == 1
        fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
        fprintf('>>> Mat''*Mat - I = \n');
        display(diffU);
        fprintf('>>> Input matrix is not unitary to the tolerance %g \n', tol);
    end;
end;

end    % of the function