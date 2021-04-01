%% *dChristMat*
% Differentiate the Christoffel matrix with respect to the stiffness coefficients and components 
% of the wavefront normal (or the slowness vector)

%%
% *Input:*
% Cij    - [6, 6] stiffness matrix in Voigt notation
% n      - [3, 1] wavefront normal or slowness
% flag   - [1, 3] array whose values equal to 1 determine which derivative arrays 
%          (dGdC, dGdN, and/or dGdNN) are to be computed

%%
% *Output:*
% dGdC   - [3, 3, 6, 6] array of the Frechet derivatives dG(i,j)/dCij(i1,j1)
% dGdN   - [3, 3, 3] array of the Frechet derivatives dG(i,j)/dn(k)
% d2GdNN - [3, 3, 3, 3] array of the Frechet derivatives dG(i,j)/(dn(k) dn(l))

%%
% Author: Vladimir Grechka 2012 - 2014

%%
function [dGdC, dGdN, d2GdNN] = dChristMat(Cij, n, flag)
%% Settings  
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

dGdC = NaN(3,3,6,6);  dGdN = NaN(3,3,3);  d2GdNN = NaN(3,3,3,3);

%% Check whether variable 'flag' is legitimate  
if isempty(flag) == 1  ||  size(flag, 1) ~= 1  ||  size(flag, 2) ~= 3
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Incorrectly specified variable ''flag'' \n');
    display(flag);
      error('>>> STOP');   
end;

%% Stiffness derivatives of the elements of the Christoffel matrix 
% G(i,j)(n) = Cij([ik],[lj])*n(k)*n(l),
% dG(i,j,i1,j1) = [d G(i,j)]/[d Cij(i1,j1)], (i, j = 1, ... 3;  i1, j1 = 1, ..., 6),
% see 'derivatives of Christoffel matrix.nb' for derivation
n1 = n(1);   n2 = n(2);   n3 = n(3);
if flag(1,1) == 1
    dGdC(:,:,1,1) = [n1^2,    0,         0;   0,           0,     0;   0,         0,       0]; 
    dGdC(:,:,1,2) = [0,       n1*n2,     0;   n1*n2,       0,     0;   0,         0,       0]; 
    dGdC(:,:,1,3) = [0,       0,     n1*n3;   0,           0,     0;   n1*n3,     0,       0]; 
    dGdC(:,:,1,4) = [0,       n1*n3, n1*n2;   n1*n3,       0,     0;   n1*n2,     0,       0]; 
    dGdC(:,:,1,5) = [2*n1*n3, 0,      n1^2;   0,           0,     0;   n1^2,      0,       0]; 
    dGdC(:,:,1,6) = [2*n1*n2, n1^2,      0;   n1^2,        0,     0;   0,         0,       0]; 
    dGdC(:,:,2,2) = [0,       0,         0;   0,        n2^2,     0;   0,         0,       0]; 
    dGdC(:,:,2,3) = [0,       0,         0;   0,           0, n2*n3;   0,     n2*n3,       0]; 
    dGdC(:,:,2,4) = [0,       0,         0;   0,     2*n2*n3,  n2^2;   0,      n2^2,       0]; 
    dGdC(:,:,2,5) = [0,       n2*n3,     0;   n2*n3,       0, n1*n2;   0,     n1*n2,       0]; 
    dGdC(:,:,2,6) = [0,       n2^2,      0;   n2^2,  2*n1*n2,     0;   0,         0,       0]; 
    dGdC(:,:,3,3) = [0,       0,         0;   0,           0,     0;   0,         0,    n3^2]; 
    dGdC(:,:,3,4) = [0,       0,         0;   0,           0,  n3^2;   0,      n3^2, 2*n2*n3]; 
    dGdC(:,:,3,5) = [0,       0,      n3^2;   0,           0,     0;   n3^2,      0, 2*n1*n3]; 
    dGdC(:,:,3,6) = [0,       0,     n2*n3;   0,           0, n1*n3;   n2*n3, n1*n3,       0]; 
    dGdC(:,:,4,4) = [0,       0,         0;   0,        n3^2, n2*n3;   0,     n2*n3,    n2^2]; 
    dGdC(:,:,4,5) = [0,       n3^2,  n2*n3;   n3^2,        0, n1*n3;   n2*n3, n1*n3, 2*n1*n2]; 
    dGdC(:,:,4,6) = [0,       n2*n3,  n2^2;   n2*n3, 2*n1*n3, n1*n2;   n2^2,  n1*n2,       0]; 
    dGdC(:,:,5,5) = [n3^2,    0,     n1*n3;   0,           0,     0;   n1*n3,     0,    n1^2]; 
    dGdC(:,:,5,6) = [2*n2*n3, n1*n3, n1*n2;   n1*n3,       0,  n1^2;   n1*n2,  n1^2,       0]; 
    dGdC(:,:,6,6) = [n2^2,    n1*n2,     0;   n1*n2,    n1^2,     0;       0,     0,       0]; 
end;
    
%% Derivatives of the elements of the Christoffel matrix G(i,j)(n) = Cij([ik],[lj])*n(k)*n(l),
% dG(i,j,k) = [d G(i,j)]/[d n(k)], (i, j, k = 1, ..., 3),
% see 'derivatives of Christoffel matrix.nb' for derivation
if flag(1,2) == 1
    dGdN(:,:,1) = [2*(n1*Cij(1,1) + n3*Cij(1,5) + n2*Cij(1,6)), ...
                   2*n1*Cij(1,6) + n3*(Cij(1,4) + Cij(5,6)) + n2*(Cij(1,2) + Cij(6,6)), ...
                   2*n1*Cij(1,5) + n3*(Cij(1,3) + Cij(5,5)) + n2*(Cij(1,4) + Cij(5,6)); ...
                   2*n1*Cij(1,6) + n3*(Cij(1,4) + Cij(5,6)) + n2*(Cij(1,2) + Cij(6,6)), ...
                   2*(n2*Cij(2,6) + n3*Cij(4,6) + n1*Cij(6,6)), ...
                   n3*(Cij(3,6) + Cij(4,5)) + n2*(Cij(2,5) + Cij(4,6)) + 2*n1*Cij(5,6); ...
                   2*n1*Cij(1,5) + n3*(Cij(1,3) + Cij(5,5)) + n2*(Cij(1,4) + Cij(5,6)), ...
                   n3*(Cij(3,6) + Cij(4,5)) + n2*(Cij(2,5) + Cij(4,6)) + 2*n1*Cij(5,6), ...
                   2*(n3*Cij(3,5) + n2*Cij(4,5) + n1*Cij(5,5))];
    dGdN(:,:,2) = [2*(n1*Cij(1,6) + n3*Cij(5,6) + n2*Cij(6,6)), ...
                   2*n2*Cij(2,6) + n3*(Cij(2,5) + Cij(4,6)) + n1*(Cij(1,2) + Cij(6,6)), ...
                   n3*(Cij(3,6) + Cij(4,5)) + 2*n2*Cij(4,6) + n1*(Cij(1,4) + Cij(5,6)); ...
                   2*n2*Cij(2,6) + n3*(Cij(2,5) + Cij(4,6)) + n1*(Cij(1,2) + Cij(6,6)), ...
                   2*(n2*Cij(2,2) + n3*Cij(2,4) + n1*Cij(2,6)), ...
                   2*n2*Cij(2,4) + n3*(Cij(2,3) + Cij(4,4)) + n1*(Cij(2,5) + Cij(4,6)); ...
                   n3*(Cij(3,6) + Cij(4,5)) + 2*n2*Cij(4,6) + n1*(Cij(1,4) + Cij(5,6)), ...
                   2*n2*Cij(2,4) + n3*(Cij(2,3) + Cij(4,4)) + n1*(Cij(2,5) + Cij(4,6)), ...
                   2*(n3*Cij(3,4) + n2*Cij(4,4) + n1*Cij(4,5))];
    dGdN(:,:,3) = [2*(n1*Cij(1,5) + n3*Cij(5,5) + n2*Cij(5,6)), ...
                   2*n3*Cij(4,5) + n2*(Cij(2,5) + Cij(4,6)) + n1*(Cij(1,4) + Cij(5,6)), ...
                   2*n3*Cij(3,5) + n2*(Cij(3,6) + Cij(4,5)) + n1*(Cij(1,3) + Cij(5,5)); ...
                   2*n3*Cij(4,5) + n2*(Cij(2,5) + Cij(4,6)) + n1*(Cij(1,4) + Cij(5,6)), ...
                   2*(n2*Cij(2,4) + n3*Cij(4,4) + n1*Cij(4,6)), ...
                   2*n3*Cij(3,4) + n2*(Cij(2,3) + Cij(4,4)) + n1*(Cij(3,6) + Cij(4,5)); ...
                   2*n3*Cij(3,5) + n2*(Cij(3,6) + Cij(4,5)) + n1*(Cij(1,3) + Cij(5,5)), ...
                   2*n3*Cij(3,4) + n2*(Cij(2,3) + Cij(4,4)) + n1*(Cij(3,6) + Cij(4,5)), ...
                   2*(n3*Cij(3,3) + n2*Cij(3,4) + n1*Cij(3,5))];
end;

%% Derivatives of the elements of the Christoffel matrix G(i,j)(n) = Cij([ik],[lj])*n(k)*n(l),
% d2G(i,j,k,l) = [d^2 G(i,j)]/[d n(k) d n(l)], (i, j, k, l = 1, ..., 3),
% see 'derivatives of Christoffel matrix.nb' for derivation
if flag(1,3) == 1
    d2GdNN(:,:,1,1) = 2*[Cij(1,1), Cij(1,6), Cij(1,5); ...
                         Cij(1,6), Cij(6,6), Cij(5,6); ...
                         Cij(1,5), Cij(5,6), Cij(5,5)];
    d2GdNN(:,:,1,2) = [2*Cij(1,6), Cij(1,2) + Cij(6,6), Cij(1,4) + Cij(5,6); ...
                       Cij(1,2) + Cij(6,6), 2*Cij(2,6), Cij(2,5) + Cij(4,6); ... 
                       Cij(1,4) + Cij(5,6), Cij(2,5) + Cij(4,6), 2*Cij(4,5)];
    d2GdNN(:,:,1,3) = [2*Cij(1,5), Cij(1,4) + Cij(5,6), Cij(1,3) + Cij(5,5); ... 
                       Cij(1,4) + Cij(5,6), 2*Cij(4,6), Cij(3,6) + Cij(4,5); ...
                       Cij(1,3) + Cij(5,5), Cij(3,6) + Cij(4,5), 2*Cij(3,5)];
    d2GdNN(:,:,2,1) = d2GdNN(:,:,1,2);               
    d2GdNN(:,:,2,2) = 2*[Cij(6,6), Cij(2,6), Cij(4,6); ... 
                         Cij(2,6), Cij(2,2), Cij(2,4); ...
                         Cij(4,6), Cij(2,4), Cij(4,4)];
    d2GdNN(:,:,2,3) = [2*Cij(5,6), Cij(2,5) + Cij(4,6), Cij(3,6) + Cij(4,5); ...
                       Cij(2,5) + Cij(4,6), 2*Cij(2,4), Cij(2,3) + Cij(4,4); ...
                       Cij(3,6) + Cij(4,5), Cij(2,3) + Cij(4,4), 2*Cij(3,4)];
    d2GdNN(:,:,3,1) = d2GdNN(:,:,1,3);               
    d2GdNN(:,:,3,2) = d2GdNN(:,:,2,3);               
    d2GdNN(:,:,3,3) = 2*[Cij(5,5), Cij(4,5), Cij(3,5); ...
                         Cij(4,5), Cij(4,4), Cij(3,4); ...
                         Cij(3,5), Cij(3,4), Cij(3,3)];
end;

end    % of the function
