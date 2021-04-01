%% *dTdR*
% Derivative of time in a homogeneous anisotropic layer with respect to the ray coordinates 

%%
% *Input:*

% r       - [3, 1] ray that generally has non-unit length   
% g       - [scalar] group velocity
% dV      - [3, 1] array of derivatives of the phase velocity calculated in 'rayAttributes'
% dg      - [2, 1] array of derivatives of the group velocity calculated in 'rayAttributes'
% flagWAA - flag indicating whether the calculations should be performed for weak
%           (flagWAA = 'WA') or strong (flagWAA = 'SA') anisotropy  

%%
% *Output:*

% dt      - [3, 1] time derivative

%%
% *Author:* Vladimir Grechka 1989 2012 

%%
function [dt] = dTdR(r, g, dV, dg, flagWAA)
%% Settings  
[~, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

rayLength = norm(r);   rayNrm = r/rayLength;

%% Compute dt/dr
if strcmp(flagWAA, 'WA') == 1
    dt = rayNrm/g - dV/g^2;
else
    dr = dVecdAngle(rayNrm);                                 % dr = dr/dtheta
    dgdr = linsolve([rayLength*dr'; rayNrm'], [dg; 0]);      % solution to eqs 1.28 - 1.30 in 
    dt = rayNrm/g - dgdr*(rayLength/g^2);                    % Obolentseva and Grechka (1989) 
end;

end    % of the function