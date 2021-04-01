%% *isTI*
% Check whether input stiffness matrix Cij is trasversely isotropic (TI) within given tolerance

%%
% *Input:*

% Cij  - [6, 6] stiffness matrix
% tol  - [scalar] tolerance of proximity to transverse isotropy  

%%
% *Output:*

% flag - [scalar] equal to 1 if Cij is TI within given tolerance and to 0 otherwise
% R    - [3, 3] rotation matrix that transforms Cij into its VTI form

%%
% *Author:* Vladimir Grechka 2012 - 2014

%%
function [flag, R] = isTI(Cij, tol)
%% Settings and defaults

% [thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
% narginTrue = nargin(thisFileName);   
% narginchk(narginTrue, narginTrue);

if nargin < 2;   tol = 1.e-6;   end                % default for tol

%% Construct two-index contractions of the stiffness tensor
Dij = zeros(3, 3);  Vij = zeros(3, 3);
[~, ~, ~, t2m] = indexesVMT(4);                     % get the indexes
for i = 1:3   
    for j = i:3
        % Dilatational stiffness tensor (Helbig, 1994, p 111; Cowin, 1989)
        Dij(i,j) = Cij(t2m(i,j), t2m(1,1)) + Cij(t2m(i,j), t2m(2,2)) + Cij(t2m(i,j), t2m(3,3));
        % Voigt stiffness tensor (Helbig, 1994, p 111; Cowin, 1989)
        Vij(i,j) = Cij(t2m(i,1), t2m(j,1)) + Cij(t2m(i,2), t2m(j,2)) + Cij(t2m(i,3), t2m(j,3));
    end   
end
[vec, val] = eig(symMat(Dij));                      % eigenvectors and eigenvalues of Dij
dval = [abs(val(2,2) - val(3,3)), abs(val(1,1) - val(3,3)), abs(val(1,1) - val(2,2))];

if sum(dval) < tol
    % If all eigenvalues of Dij are equal, its eigenvectors are indeterminable 
    % -> Use the Voigt tensor instead
    [vec, val] = eig(symMat(Vij));                   % eigenvectors and eigenvalues of Vij
    dval = [abs(val(2,2) - val(3,3)), abs(val(1,1) - val(3,3)), abs(val(1,1) - val(2,2))];
end

%% Get direction of the symmetry axis
% Since two eigenvalues of Dij are equal in TI media, a numerically zero difference between 
% the eigenvalues is a necessary condition for TI
[valsort, isort] = sort(dval, 'descend');

if valsort(3) < tol
    % If medium is TI, the third column of R points in direction of the symmetry axis     
    R = vec(:,isort);   
    CijRot = bond(Cij, R);                          % rotation to VTI if 'Cij' is TTI
    
    % Check whether 'CijRot' is, indeed, a VTI matrix
    check1 = sum(abs(CijRot(1:3, 4))) + ...         % off-diagonal zeros
             sum(abs(CijRot(1:4, 5))) + ...
             sum(abs(CijRot(1:5, 6)));   
    check2 = abs(CijRot(1, 1) - CijRot(2, 2)) + ... % constraints on the diagonal elements
             abs(CijRot(1, 3) - CijRot(2, 3)) + ...
             abs(CijRot(4, 4) - CijRot(5, 5)) + ...
             abs(CijRot(1, 2) + 2*CijRot(6, 6) - CijRot(1, 1));
    if check1 + check2 < tol                                     
        flag = 1;                                   % 'Cij' has TI symmetry
    else
        flag = 0;                                   % The symmetry of 'Cij' is lower than TI
        R = eye(3, 3);
    end
else                                                % The symmetry of 'Cij' is lower than TI
    flag = 0;                                       
    R = eye(3, 3);
end

end    % of the function