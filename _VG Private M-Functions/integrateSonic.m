%% *integrateSonic*
% Integrate sonic log 

%%
% *Input:*

% depth           - [:, 1] array of sonic depths
% sonicDepthUnits - [string] units of the sonic depth
% zRec            - [:, 1] array of receiver depths  
% recLengthUnits  - [string] units of the receiver depth
% vel             - [:, 1] array of sonic velocities
% velLengthUnits  - [string] length units of the velocities 
% velTimeUnits    - [string] time units of the velocities 
% izSou           - [scalar] the number of receiver turned into a source
% updownFlag      - [scalar] flag equal to +1 or -1 to indicate up- or down-propagating wave
% outTimeUnits    - [string] time units of the output

%%
% *Output:*

% timeSonic       - [size(zRec, 1), 1] array of times of the integrated sonic log 

%%
% *Author:* Vladimir Grechka 2012 2013    

%%
function [timeSonic] = integrateSonic( ...
    depth, sonicDepthUnits, zRec, recLengthUnits, vel, velLengthUnits, velTimeUnits, ...
    izSou, updownFlag, outTimeUnits)
%% Settings 
[~, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

%% Change the length units to [m] and the time units - to [outTimeUnits]
depth = u2u(depth, sonicDepthUnits, 'm', []);
zRec  = u2u(zRec,  recLengthUnits,  'm', []);
vel   = u2u(vel,   velLengthUnits,  'm', []);
vel   = u2u(vel,   outTimeUnits,    velTimeUnits, []);

%% Replace the missing (NaN) velocities with the mean velocity 
vMean = mean(vel(~isnan(vel)));
vel(isnan(vel)) = vMean;  

%% Fill possibly missing beginning or end of the log
dz = depth(2,1) - depth(1,1);
if zRec(1) < depth(1,1);     
    zz1  = zRec(1) : dz : depth(1,1);
    vel1 = vMean*ones(length(depth1), 1); 
else
    zz1 = [];   vel1 = [];   
end;

if depth(end,1) < zRec(end);
    zz2 = depth(end,1) : dz : zRec(end);
    vel2 = vMean*ones(length(depth2), 1); 
else
    zz2 = [];   vel2 = [];   
end;

zzFull  = cat(1,  zz1, depth,  zz2);
velFull = cat(1, vel1,   vel, vel2);

%% Map the receiver depths onto the sonic depths
[~, indBeg] = min(abs(zRec(1)   - zzFull));
[~, indEnd] = min(abs(zRec(end) - zzFull));

%% Integrate the sonic
time = zeros(indEnd - indBeg, 1);
for i = 2:(indEnd - indBeg)
    time(i,1) = time(i-1,1) + dz/velFull(indBeg + i - 1);
end;
zarray = zRec(1) : dz : zRec(end);
timeSonic = interp1(zarray(1:min([length(zarray), length(time)])), ...
                      time(1:min([length(zarray), length(time)])), zRec', 'linear', 'extrap');
timeSonic = updownFlag*(timeSonic - timeSonic(izSou))';

end    % of the function