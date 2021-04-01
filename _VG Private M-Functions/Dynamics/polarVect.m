%% *polarVect*
% Solve the Christoffel equation for the polarization vector 

%%
% *Input:*

% Cij - [6, 6] stiffness matrix
% p   - [3, 1] slowness vectors

%%
% *Output:*

% U   - [3, 1] polarization vector (column)
%
%%
% *Author:* Geser Dugarov 2018

%%
function [U] = polarVect(Cij, p)

G = christoffel(Cij, p, p);
[U, V] = eig(G);
temp = diag(V);
temp = abs(temp - 1);
U = U(:,temp < 1.0e-10);

end    % of the function

%%
