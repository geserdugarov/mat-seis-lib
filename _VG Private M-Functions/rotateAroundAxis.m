%% *rotateAroundAxis*
% Construct matrix R describing rotation around an axis  

%% 
% *Input:* 

% rotAxis    - [3, 1] vector describing the axis of rotation
% rotAngle   - [scalar] rotation angle (in deg or rad)
% angleUnits - [string] 'deg' or 'rad'

%%
% *Output:*

% Rot      - [3, 3] rotation matrix

%%
% *Author:* Vladimir Grechka 2012 - 2014

function [Rot] = rotateAroundAxis(rotAxis, rotAngle, angleUnits)
%% Settings and defaults
[~, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

rotAxis  = rotAxis/norm(rotAxis);               % normalize 'rotAxis'
rotAngle = u2u(rotAngle, angleUnits, 'rad');    % make 'rotAngle' in rad

%% Rotation matrix 
%  (Korn and Korn, 1984, p 448, eq 14.10-6)
Rot = cos(rotAngle)*eye(3) + ...
    (1 - cos(rotAngle))* ...
        [rotAxis(1)*rotAxis(1), rotAxis(1)*rotAxis(2), rotAxis(1)*rotAxis(3); ...
         rotAxis(2)*rotAxis(1), rotAxis(2)*rotAxis(2), rotAxis(2)*rotAxis(3); ...
         rotAxis(3)*rotAxis(1), rotAxis(3)*rotAxis(2), rotAxis(3)*rotAxis(3)] + ...
    sin(rotAngle)* ...
        [    0,      -rotAxis(3),  rotAxis(2); ...
          rotAxis(3),    0,       -rotAxis(1); ...
         -rotAxis(2), rotAxis(1),     0];
             
end  % of the function
