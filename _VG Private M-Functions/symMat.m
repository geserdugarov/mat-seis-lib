%% *symMat*
% Symmetrize a square matrix by copying its upper-right corner into the low-left corner 

%%
% *Input:*

% Minp - square matrix whose relevant elements are located on the main diagonal and above it 

%%
% *Output:*

% Msym - symmetric square matrix generated from Minp

%%
% *Author:* Vladimir Grechka 2012

%%
function [Msym] = symMat(Minp); 
%% Settings and checks
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

% Check whether Minp is square
[n, m] = size(Minp);
if n ~= m
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    display(Minp)
    fprintf('>>> Input matrix ''Minp'' is not square \n \n')
    error('>>> STOP');
end; 

%% Fill the symmetric part
Msym = Minp;
for i = 1:n
    for j = 1:i-1
        Msym(i,j) = Msym(j,i);
    end;
end;

end    % of the function