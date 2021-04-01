%==================================================================================================
% Set planar model interfaces 
%
% function [fault, noFlt] = setFault;
%
% Input: none
%
% Output: 
%   fault - [5, noFlt] fault array  
%   noInt - the number of faults
%       Array 'fault' contains columns with elements
%       [xyzFlt(1), xyzFlt(2), xyzFlt(3), nrmFlt(1), nrmFlt(2)],
%       in which xyzFlt(i) are the coordinates of a point at the fault and nrmFlt(i) are the 
%       horizontal components of a unit normal to the fault
%           Faults are described by equations dot(nrmFlt, x - xyzFlt) = 0 with 
%               nrmFlt(3) = sqrt( 1 - (nrmInt(1)^2 + nrmInt(2)^2) ) > 0
%           Faults can intersect interfaces and each other and can be inputed in any order
%           Empty 'depth' array means no faults 
%
% Copyright (C) Vladimir Grechka 2012
%==================================================================================================

function [fault, noFlt] = setFault;

%==================================================================================================
functionInfo = functions(@setFault);   

% Fault depths, dips, and azimuths
%depth = [0.15];   % four-layer model
depth = [];   
noFlt = length(depth);                              % the number of faults
dip = u2u(80.0*ones(1,noFlt), 'deg', 'rad', 0);         % fault dip
azm = u2u(180.0*ones(1,noFlt), 'deg', 'rad', 0);         % fault azimuth 

if isempty(depth) == 1;   fault = zeros(5,0);   return;   end;

% Set the output 'fault' array
fault = zeros(5,noFlt);                             % intial setting
fault(1:2,noFlt) = [0.3, 0];                               % the coordinate origin for the faults
fault(  3,:) = depth;                               % fault depths
fault(  4,:) = sin(dip).*cos(azm);                  % fault normal component nrmFlt(1)
fault(  5,:) = sin(dip).*sin(azm);                  % fault normal component nrmFlt(2)

%==================================================================================================
%==================================================================================================

