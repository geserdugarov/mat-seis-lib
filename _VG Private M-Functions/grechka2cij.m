%% *grechka2cij*
% Calculate stiffness matrix in a monoclinic medium from its Grechka's (2000) parameters 

%%
% *Input:*

% Ani - [1, 12] vector of Grechka's (2000) anisotropy parameters arranged as
%       [Vp0, Vs0, epsilon1, epsilon2, delta1, delta2, delta3, gamma1, gamma2, zeta1, zeta2, zeta3] 

%%
% *Output:*

% Cij - [6, 6] monoclinic stiffness matrix

%%
% *Author:* Vladimir Grechka 2000 2012

%%
function [Cij] = grechka2cij(Ani)
%% Settings 
[~, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

%% Compute the anisotropy coefficients
Cij = tsvankin2cij(Ani(1:9));

Cij(3,6) =   Ani(12)*Cij(3,3);              Cij(6,3) = Cij(3,6);
Cij(1,6) = 2*Ani(10)*Cij(3,3) + Cij(3,6);   Cij(6,1) = Cij(1,6);      
Cij(2,6) = 2*Ani(11)*Cij(3,3) + Cij(3,6);   Cij(6,2) = Cij(2,6);   

end    % of the function