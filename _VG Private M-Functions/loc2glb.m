%% *loc2glb*
% Construct rotation matrix from local coordinate frame (x1, y1, z1) to global frame (x, y, z)

%%
% *Input:*

% angles - [3, 1] array of three angles (in rad) specifying rotation:
%          . azim - azimuth of the local z1-axis with respect to the global x-axis
%          . tilt - tilt of the local z1-axis with respect to the global z-axis
%          . azx1 - azimuth of the local x1-axis with respect to the global x-axis

%%  
% *Output:*

% Rot    - [3, 3] rotation matrix

%%
% *Author:* Vladimir Grechka 1989 1998 2012

%%
function [Rot] = loc2glb(angles)
%% Settings 
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

tol = 1.e-8;                % default for the tolerance

%% Construct the rotation matrix Rot
azim = angles(1,1);   tilt = angles(2,1);   azx1 = angles(3,1);
sa = sin(azim);   ca = cos(azim);   sb = sin(tilt);   cb = cos(tilt);
sd = sin(azx1);   cd = cos(azx1);   cda = ca*cd + sa*sd;

if abs(sb) < tol  ||  abs(cda) < tol
    ct = 0;   st = 1;       % avoid division by zero
else
    tt = -cb/(sb*cda);   ct = -1/sqrt(1 + tt^2);   st = tt*ct;
end;

x1 = [cd*st, sd*st, ct];   x3 = [ca*sb, sa*sb, cb];   x2 = cross(x3, x1);
Rot = [x1; x2; x3];

%% Make sure that Rot is unitary
if isUnitary(Rot) == 0      % make sure that Rot is unitary
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
      error('>>> STOP');
end;

end    % of the function