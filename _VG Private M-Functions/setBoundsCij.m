%% *setBoundsCij*
% Set the lower and upper bounds on the components of elastic stiffness tensor 

%%
% *Input:*  

% Cij       - [6, 6] elastic stiffness matrix
% indexCij  - [1, :] vector of indexes of the stiffness elements whose bounds are to be computed
% boundsCij - [1, 5] vector of bounds (in fractions of unity) whose components correspond to
%             . boundsCij(1) - the (lambda + 2*mu)-type elements Cij(i,i), (i = 1, 2, 3)
%             . boundsCij(2) - the lambda-type elements Cij(i,j), (i, j = 1, 2, 3)
%             . boundsCij(3) - the mu-type elements Cij(i,i), (i = 4, 5, 6)
%             . boundsCij(4) - the upper-right corner elements Cij(i,j), (i = 1, 2, 3; j = 4, 5, 6)
%             . boundsCij(5) - the lower-right corner elements Cij(i,j), (i, j = 4, 5, 6;  i ~= j)

%%
% *Output:*

% CijL      - [1, length(indexCij)] the lower bound of Cij 
% CijU      - [1, length(indexCij)] the upper bound of Cij 

%%
% *Author:* Vladimir Grechka 2012

%%
function [CijL, CijU] = setBoundsCij(Cij, indexCij, boundsCij)
%% Settings
[~, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

noCij = length(indexCij);

meanCij = abs([mean([Cij(1,1), Cij(2,2), Cij(3,3)]), ...    % get the mean values of 
               mean([Cij(1,2), Cij(1,3), Cij(2,3)]), ...    % the stiffnesses 
               mean([Cij(4,4), Cij(5,5), Cij(6,6)])]);
meanCij(4:5) = meanCij(2:3);

%% Calculate the bounding stiffness elements
CijL = NaN(1, noCij);    CijU = NaN(1, noCij);
for i = 1:noCij
    % Reminder of the v2m conversion rule in function 'indexesVMT'
    %   [1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 5, 5, 6; ...
    %    1, 2, 3, 4, 5, 6, 2, 3, 4, 5, 6, 3, 4, 5, 6, 4, 5, 6, 5, 6, 6];
    %    -------------------------------------------------------------
    %    1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21
    if     indexCij(i) ==  1 || indexCij(i) ==  7 || indexCij(i) == 12 
        ij = 1;     % the (lambda + 2*mu)-type elements Cij(i,i), (i = 1, 2, 3)
    elseif indexCij(i) ==  2 || indexCij(i) ==  3 || indexCij(i) == 8 
        ij = 2;     % the lambda-type elements Cij(i,j), (i, j = 1, 2, 3)
    elseif indexCij(i) == 16 || indexCij(i) == 19 || indexCij(i) == 21 
        ij = 3;     % the mu-type elements Cij(i,i), (i = 4, 5, 6)
    elseif indexCij(i) == 17 || indexCij(i) == 18 || indexCij(i) == 20 
        ij = 5;     % the lower-right corner elements Cij(i,j), (i, j = 4, 5, 6;  i ~= j)
    else
        ij = 4;     % the upper-right corner elements Cij(i,j), (i = 1, 2, 3;  j = 4, 5, 6)
    end;

    %% Compute the bounds
    CijL(1,i) = (-1)*boundsCij(ij)*meanCij(ij);
    CijU(1,i) =      boundsCij(ij)*meanCij(ij);    
end;

end    % of the function