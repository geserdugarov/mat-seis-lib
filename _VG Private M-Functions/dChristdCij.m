%% *dChristdCij*
% Differentiate determinant of the Christoffel matrix with respect to stiffness coefficients 
% (eq A-11 in Grechka and Duchkov, 2011)

%%
% *Input:*

% Cij - [6, 6] stiffness matrix
% p   - [3, 1] wavefront normal or slowness vector  
% V   - [scalar] phase velocity 

%%
% *Output:*

% dF  - [6, 6] matrix of the Frechet derivatives dF/dCij(i1,j1), where
%       F = det[G(i,j) - delta(i,j)],      G(i,j)(p) = Cij([ik],[lj])*p(k)*p(l) or
%       F = det[G(i,j) - V^2^delta(i,j)],  G(i,j)(n) = Cij([ik],[lj])*n(k)*n(l)

%%
% *Author:* Vladimir Grechka 2012 - 2014

%%
function [dF] = dChristdCij(Cij, p, V)
%% Settings  
[~, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

%% Construct the Christoffel matrix G
% G is applicable to either the slowness and wavefront normal vectors
%   if p is the slowness vector, V*norm(p) = 1
%   if p is the wavefront normal, that is, p = n, V*norm(p) = V*norm(n) = V
G = christoffel(Cij, p, p) - (V*norm(p))^2*eye(3);

%% Compute derivatives of the elements of the Christoffel matrix 
% dG(i,j, i1,j1) = [d G(i,j)]/[d Cij(i1,j1)], (i, j = 1, ... 3;  i1, j1 = 1, ..., 6)
[dG, ~, ~] = dChristMat(Cij, p, [1, 0, 0]); 

%% Calculate the derivatives dF/dCij
dF = zeros(6, 6);
for i1=1:6;   
    for j1=i1:6
        dF(i1,j1) = dDetdX(G, dG(:,:,i1,j1));
        
        if j1 > i1
            dF(j1,i1) = dF(i1,j1);    % symmetrize dF 
        end;
    end;   
end;

end    % of the function