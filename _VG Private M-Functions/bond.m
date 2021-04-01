%% *bond*
% Compute the Bond-transformation matrix and rotate a stiffness matrix with it

%%
% *Input:*

% Cij  - [6, 6] stiffness matrix
% Rot  - [3, 3] rotation matrix (supposedly unitary)
% flag - [scalar] indicator of whether only the Bond matrix is to be computed (flag = 0) or  
%        both the Bond matrix and the Bond transformation of Cij (flag is absent or flag ~= 0)  

%%
% *Output:*

% Crot - [6, 6] rotated stiffness matrix
% B    - [6, 6] the Bond matrix

%%
% *Author:* Vladimir Grechka 1998 2012

%%
function [Crot, B] = bond(Cij, Rot, flag)
%% Settings 
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
% narginTrue = nargin(thisFileName);   
% narginchk(narginTrue, narginTrue);

% Check whether the rotation matrix is unitary
if isUnitary(Rot)  == 0   
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Continue computations with non-unitary rotation matrix? -- PAUSE \n');
    pause;
end;

%% Construct the Bond matrix B from the rotation matrix Rot (Auld, 1990, eq 3.34)
B = [Rot(1,1)^2, Rot(1,2)^2, Rot(1,3)^2, ...
        Rot(1,2)*Rot(1,3), Rot(1,3)*Rot(1,1), Rot(1,1)*Rot(1,2); ...
     Rot(2,1)^2, Rot(2,2)^2, Rot(2,3)^2, ...
        Rot(2,2)*Rot(2,3), Rot(2,3)*Rot(2,1), Rot(2,1)*Rot(2,2); ...
     Rot(3,1)^2, Rot(3,2)^2, Rot(3,3)^2, ...
        Rot(3,2)*Rot(3,3), Rot(3,3)*Rot(3,1), Rot(3,1)*Rot(3,2); ...
     2*Rot(2,1)*Rot(3,1), 2*Rot(2,2)*Rot(3,2), 2*Rot(2,3)*Rot(3,3), ...
       Rot(2,2)*Rot(3,3) + Rot(2,3)*Rot(3,2), Rot(2,1)*Rot(3,3) + Rot(2,3)*Rot(3,1), ...
       Rot(2,2)*Rot(3,1) + Rot(2,1)*Rot(3,2); ...
     2*Rot(3,1)*Rot(1,1), 2*Rot(3,2)*Rot(1,2), 2*Rot(3,3)*Rot(1,3), ...
       Rot(1,2)*Rot(3,3) + Rot(1,3)*Rot(3,2), Rot(1,3)*Rot(3,1) + Rot(1,1)*Rot(3,3), ...
       Rot(1,1)*Rot(3,2) + Rot(1,2)*Rot(3,1); ...
     2*Rot(1,1)*Rot(2,1), 2*Rot(1,2)*Rot(2,2), 2*Rot(1,3)*Rot(2,3), ...
       Rot(1,2)*Rot(2,3) + Rot(1,3)*Rot(2,2), Rot(1,3)*Rot(2,1) + Rot(1,1)*Rot(2,3), ...
       Rot(1,1)*Rot(2,2) + Rot(1,2)*Rot(2,1)];

if nargin == 3  &&  flag == 0  
    Crot = NaN(6, 6);
else
    % The Bond transformation
    Crot = B'*Cij*B; 
end;

end    % of the function