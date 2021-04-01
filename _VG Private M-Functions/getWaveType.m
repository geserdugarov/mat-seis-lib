%% *getWaveType*
% Compute the wave type in the current layer relying on the polarization continuity

%%
% *Input:*

% Uprev    - [3, 1] polarization vector in the previous layer
% Ucurr    - [3, 3] polarization matrix in the current layer

%%
% *Output:*

% waveType - [scalar] equal to 1 (for the P-wave), or to 2 (for the S1-wave), 
%            or to 3 (for the S2-wave)

%%
% *Author:* Vladimir Grechka 2013 

%%
function [waveType] = getWaveType(Uprev, Ucurr)
%% Settings 
[~, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

%% Get the dot products between polarization vectors in the current layer and the polarization vectors in the 
tmp = abs(dot(Ucurr, [Uprev, Uprev, Uprev], 1))

%% Calculate the current 'waveType'
waveType = find(tmp == max(tmp));

end    % of the function