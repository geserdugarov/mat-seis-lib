%% *rayTracingOptim*
% General optimization procedure for two-point ray tracing in layered anisotropic media

%%
% *Input:*

% xy      - [2*(noSeg-1), 1] array of the pairs of x- and y-coordinates of intersections of a 
%           given ray with the model interfaces and faults
% rayType - [noSeg, 1] array of ray types (1 for P, 2 for S1 or SV, and 3 for S2 or SH wave) 
%           along each ray segment
%   noSeg - the number of ray segments
% xSou    - [3, 1] source coordinates
% xRec    - [3, 1] receiver coordinates 
% zInt    - [3, noSeg-1] array of the coefficients of interfaces given by 
%           z = x*zInt(1) + y*zInt(2) + zInt(3)
%           (*) The interface depth at the i-th ray intresection is 
%               z = dot(zInt(1:2,i), xy(2*i-1 : 2*i, 1)) + zInt(3,i)
% Cij     - [6, 6, noSeg] stiffness matrices of layers arranged along the ray trajectory
% flagWAA - flag indicating whether calculations should be performed for weak
%           (flagWAA = 'WA') or strong (flagWAA = 'SA') anisotropy  

%%
% *Output:*

% time    - time along a ray
% dtdx    - [2*(noSeg-1), 1] gradient of time [dtime/dx(i); dtime/dy(i)] for i = [1:noSeg-1]  
% d2tdxy  - [2*(noSeg-1), 2*(noSeg-1)] Hessian matrix

%%
% *Author:* Vladimir Grechka 2012 - 2014

%%
function [time, dtdx, d2tdxy] = rayTracingOptim(xy, rayType, xSou, xRec, zInt, Cij, flagWAA)
%function [time, dtdx] = rayTracingOptim(xy, rayType, xSou, xRec, zInt, Cij, flagWAA)
%% Settings 
% [~, thisFileName, ~] = fileparts(mfilename('fullpath'));
% narginTrue = nargin(thisFileName);   
% narginchk(narginTrue, narginTrue);
noSeg = size(rayType, 1);

%% Get the ray attributes
[xyz, tOut, ~, ~, ~, ~, ~, dVOut, gOut, dgOut, ~, ~, flagOut] = ...
    rayAttributes(xy, rayType, xSou, xRec, zInt, Cij, flagWAA);
if flagOut == 0
    time = 1.e+10;                      % The ray attributes have not been found -- 
    dtdx = zeros(2*(noSeg-1), 1);       % force a quick exit
    if nargout > 2
        d2tdxy = eye(2*(noSeg-1));
    end
    return;
end

%% Time and its gradient and Hessian for each ray segment
time = 0;   dtdxSeg = zeros(3, noSeg);   d2tdxySeg = zeros(3, 3, noSeg);
for iseg = 1:noSeg
    raySeg = (xyz(iseg+1,:) - xyz(iseg,:));    rayLength = norm(raySeg);   
    time = time + tOut(1,iseg);
    dtdxSeg(:,iseg) = dTdR(raySeg', gOut(iseg), dVOut(:,iseg), dgOut(:,iseg), flagWAA);
    if nargout > 2
        rayNrm = raySeg/rayLength;
        d2tdxySeg(:,:,iseg) = (eye(3) - rayNrm'*rayNrm)/(rayLength*gOut(iseg));
    end
end

%% Gradient and Hessian of the time along the entire ray trajectory
dtdx = zeros(2*(noSeg-1), 1);
if nargout > 2
    d2tdxy = zeros(2*(noSeg-1), 2*(noSeg-1));
end
for iseg = 1:(noSeg-1)
    mat1 = [eye(2), zInt(1:2,iseg)];
    dtdx(2*iseg-1 : 2*iseg, 1) = mat1*(dtdxSeg(:,iseg) - dtdxSeg(:,iseg+1));
    if nargout > 2
        d2tdxy(2*iseg-1 : 2*iseg, 2*iseg-1 : 2*iseg) = ...
                             mat1*(d2tdxySeg(:,:,iseg) + d2tdxySeg(:,:,iseg+1))*mat1';
    end
end

if nargout > 2
    for iseg = 1:(noSeg-2)
        mat2 = -[eye(2), zInt(1:2,iseg)]*d2tdxySeg(:,:,iseg+1)*[eye(2), zInt(1:2,iseg+1)]';
        d2tdxy(2*iseg-1 : 2*iseg,   2*iseg+1 : 2*iseg+2) = mat2;
        d2tdxy(2*iseg+1 : 2*iseg+2, 2*iseg-1 : 2*iseg)   = mat2';
    end
end

end    % of the function
