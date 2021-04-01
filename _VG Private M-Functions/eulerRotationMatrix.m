%% *eulerRotationMatrix*
% Construct Euler rotation matrix  

%% 
% *Input:* 

% angles     - [3, 1] vector of Euler angles alpha, beta, gamma
% angleUnits - [string] 'deg' or 'rad'

%%
% *Output:*

% Rot      - [3, 3] rotation matrix

%%
% *Author:* Vladimir Grechka 2014

function [Rot] = eulerRotationMatrix(angles, angleUnits)
%% Settings and defaults
[~, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

ang = u2u(angles, angleUnits, 'rad');    % make 'angles' in rad

%% Euler rotation matrix 
%  (http://en.wikipedia.org/wiki/Euler_angles)
s = sin(ang);
c = cos(ang);
Rot = [c(2),      -c(3)*s(2),                   s(2)*s(3); ...
       c(1)*s(2),  c(1)*c(2)*c(3) - s(1)*s(3), -c(3)*s(1) - c(1)*c(2)*s(3); ...
       s(1)*s(2),  c(1)*s(3) + c(2)*c(3)*s(1),  c(1)*c(3) - c(2)*s(1)*s(3)];
                
end  % of the function

%%