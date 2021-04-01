%% *setInterface*
% Specify planar model interfaces 

%% 
% *NB:*
%
% * This script needs to be tailored to each specific data set

%%
% *Input:* none

%%
% *Output:*

% interface - [5, noInt] interface array  
% noInt     - the number of interfaces
%             (*) Array 'interface' contains columns with elements
%                 [xyzInt(1), xyzInt(2), xyzInt(3), nrmInt(1), nrmInt(2)],
%                 in which xyzInt(i) are the coordinates of a point at the interface and 
%                 nrmInt(i) are the horizontal components of a unit normal to the interface
%             (*) Interfaces are described by equations dot(nrmInt, x - xyzInt) = 0 with 
%                 nrmInt(3) = sqrt( 1 - (nrmInt(1)^2 + nrmInt(2)^2) ) > 0
%             (*) Interfaces should not intersect in the volume covered by rays and are 
%                 ordered in the increasing depths xyzInt(3) 
%             (*) Empty 'depth' array implies a homogeneous space and noInt = 0

%%
% *Author:* Vladimir Grechka 2012

%%
function [interface, noInt] = setInterface

% Settings  
% [thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
% narginTrue = nargin(thisFileName);   
% narginchk(narginTrue, narginTrue);

%% Interface depths, dips, and azimuths
global modelflag;
if modelflag==1
    depth = [1000 1500];
elseif modelflag==2
    depth = [1000 1500];
elseif modelflag==3
    depth = [1000 1500];
elseif modelflag==4
    depth = [1300 1580];
elseif modelflag==5
    depth = [100, 740, 1010, 1700, 2230, 2820, 2960];
elseif modelflag==6
    depth = [10, 220, 480, 640, 770, 1050, 1300, 1580, 1620];
elseif modelflag==7
    depth = [400, 1000, 1700, 2220];
end
noInt = length(depth);                                 % the number of interfaces
dip = u2u(0.0*ones(1,noInt), 'deg', 'rad', 0);         % interface dip
azm = u2u(0.0*ones(1,noInt), 'deg', 'rad', 0);         % interface azimuth 
% depth = [1000, 1900, 2500];
%noInt = length(depth);                            % the number of interfaces
% dip = u2u([2,  3, 4], 'deg', 'rad', 0);          % interface dip
% azm = u2u([10, 5, -5], 'deg', 'rad', 0);         % interface azimuth

if isempty(depth) == 1;   interface = zeros(5,0);   return;   end
   
%% Set the output 'interface' array
interface = zeros(5,noInt);                         % intial setting
interface(1:2,noInt) = 0;                           % the coordinate origin for the interfaces
interface(  3,:) = depth;                           % interface depths
interface(  4,:) = sin(dip).*cos(azm);              % interface normal component nrmInt(1)
interface(  5,:) = sin(dip).*sin(azm);              % interface normal component nrmInt(2)

end    % of the function

% depthPlane = [0.,   0.5,   1.0];
% dipPlane   = [0.,  25.0,   5.0];
% azmPlane   = [0.,  10.0, -20.0];