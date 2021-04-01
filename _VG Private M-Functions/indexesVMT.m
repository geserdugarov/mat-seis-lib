%% *indexesVMT*
% Convert stiffness matrix in Voigt notation into its row-vector or tensor representations 
% or restore the Voigt matrix back

%%
% *Input:*

% flag - [scalar] equal to 1, 2, 3, or 4 to indicate which of the four outputs is to be computed  

%%
% *Output:*

% v2m  - [2, 21] array converting the [1, 21] row-vector of stiffnesses into the [6, 6] Voigt
%        stiffness matrix     
% m2v  - [6, 6] array converting the [6, 6] Voigt stiffness matrix into the [1, 21] row-vector 
%        of stiffnesses
% m2t  - [6, 2] array converting the [6, 6] Voigt stiffness matrix into the [3, 3, 3, 3]
%        stiffness tensor
% t2m  - [3, 3] array converting the [3, 3, 3, 3] stiffness tensor into the [6, 6] Voigt
%        stiffness matrix 

%%
% *Author:* Vladimir Grechka 2012

%%
function [v2m, m2v, m2t, t2m] = indexesVMT(flag)
%% Settings and defaults
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

v2m = NaN(2, 21);   m2v = NaN(6, 6);    m2t = NaN(6, 2);    t2m = NaN(3, 3);

%% Checks
if isempty(flag) == 1   
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Empty variable ''flag'' \n \n');
      error('>>> STOP')     
end;

if (isequal(flag, 1) + isequal(flag, 2) + isequal(flag, 3) + isequal(flag, 4)) ~= 1
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Illegitimate value of variable flag = %g \n', flag);
    fprintf('    It should be equal to either 1, 2, 3, or 4 \n \n');
      error('>>> STOP')     
end;

%% Create output arrays
if flag == 1
    % Row-vector-to-matrix conversion: 1 -> 11, ..., 6 -> 16, 7 -> 22, ..., 21 -> 66 
    % Cm(v2m(1,k), v2m(2,k)) = Cv(k),  (k = 1, ..., 21)
    v2m = [1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 5, 5, 6; ...
           1, 2, 3, 4, 5, 6, 2, 3, 4, 5, 6, 3, 4, 5, 6, 4, 5, 6, 5, 6, 6];
end;

if flag == 2
    % Voigt matrix-to-row-vector conversion: 11 -> 1, ..., 16 -> 6, 22 -> 7, ..., 66 -> 21 
    % Cv(m2v(i,j)) = Cm(i,j),  (i, j = 1, ..., 6)
    m2v = [1,  2,  3,  4,  5,  6; ... 
           2,  7,  8,  9, 10, 11; ... 
           3,  8, 12, 13, 14, 15; ...
           4,  9, 13, 16, 17, 18; ... 
           5, 10, 14, 17, 19, 20; ... 
           6, 11, 15, 18, 20, 21];
end;

if flag == 3
    % Voigt matrix-to-tensor index conversion: 1 -> 11, 2 -> 22, 3 -> 33, 4 -> 23, 5 -> 13, 6 -> 12
    % Ct(m2t(i,1), m2t(i,2), m2t(j,1), m2t(j,2)) = Cm(i,j),  (i, j = 1, ..., 6)
    m2t = [1, 1; ... 
           2, 2; ... 
           3, 3; ... 
           2, 3; ... 
           1, 3; ... 
           1, 2];
end;

if flag == 4
    % Tensor-to-matrix conversion: 11 -> 1, 22 -> 2, 33 -> 3, 23 -> 4, 13 -> 5, 12 -> 6 
    % Cm(t2m(i1,i2), t2m(j1,j2)) = Ct(i1,i2,j1,j2),  (i1, i2, j1, j2 = 1, 2, 3)
    t2m = [1, 6, 5; ... 
           6, 2, 4; ... 
           5, 4, 3];
end;

end    % of the function 