%% *readSonic*
% Read a spreadsheet containing sonic data 

%%
% *Input:*

% fullFileName             - [string] full name of file containing sonic data in .xls 
%                            or .xlsx format
% colDepth                 - [scalar] column containing the depth        
% inpDepthUnits            - [string] units of depths in 'colDepth'
% outDepthUnits            - [string] units of depths in output structure 'sonic'  
% colVp                    - [scalar] column containing the P-wave velocity   
% colVs1                   - [scalar] column containing the S1-wave velocity  
% colVs2                   - [scalar] column containing the S2-wave velocity  
% inpVelocityLengthUnits   - [string] length units of the velocities in columns 'colVp', 
%                            'colVs1', and 'colVs2'  
% inpVelocityTimeUnits     - [string] time units of the velocities in columns 'colVp', 
%                            'colVs1', and 'colVs2'  
% outVelocityLengthUnits   - [string] length units for the velocities in output structure 'sonic'                         
% outVelocityTimeUnits     - [string] time units for the velocities in output structure 'sonic'
% colQp                    - [scalar] column containing the P-wave slowness
% colQs1                   - [scalar] column containing the S1-wave slowness
% colQs2                   - [scalar] column containing the S2-wave slowness
% inpSlownessLengthUnits   - [string] length units of the slownesses in columns 'colQp', 
%                            'colQs1', and 'colQs2'  
% inpSlownessTimeUnits     - [string] time units of the slownesses in columns 'colQp', 
%                            'colQs1', and 'colQs2'
% missingData              - [scalar] indicating a missing value
% flagNaN                  - [scalar] flag equal to 1 to replace missing values with NaNs 
%                            and ~= 1 to leave them as they are    
% (*) Empty column numbers for the depth, velocities and slownesses result in empty fields
%     of structure 'sonic'

%%
% *Output:*

% sonic                    - structure containing the following fields:        
%     sonic.fullFileName   - the name of file containing the original sonic
%     sonic.Depth          - depth of sonic in 'sonic.DepthUnits'
%     sonic.DepthUnits     - units of 'sonic.Depth'
%     sonic.Vp             - sonic P-wave velocity 
%     sonic.Vs1            - sonic S1-wave velocity 
%     sonic.Vs2            - sonic S2-wave velocity 
%     sonic.VelLengthUnits - the length units of velocities 'sonic.Vp', 'sonic.Vs1',
%                            and 'sonic.Vs2'
%     sonic.VelTimeUnits   - the time units of velocities 'sonic.Vp', 'sonic.Vs1',
%                            and 'sonic.Vs2'

%%
% *Author:* Vladimir Grechka 2012 - 2014    

%%
function [sonic] = ...
    readSonic(fullFileName, colDepth, inpDepthUnits, outDepthUnits, ...
              colVp, colVs1, colVs2, inpVelocityLengthUnits, inpVelocityTimeUnits, ...
                                     outVelocityLengthUnits, outVelocityTimeUnits, ...
              colQp, colQs1, colQs2, inpSlownessLengthUnits, inpSlownessTimeUnits, ...
              missingData, flagNaN)
%% Settings
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

sonic.fullFileName   = fullFileName;
sonic.DepthUnits     = outDepthUnits;   
sonic.VelLengthUnits = outVelocityLengthUnits;   
sonic.VelTimeUnits   = outVelocityTimeUnits;

sonic.Depth = [];   sonic.Vp  = [];   sonic.Vs1 = [];   sonic.Vs2 = []; 
tol = 1.e-12;

%% Read the sonic
[num, ~, ~] = xlsread(sonic.fullFileName);

if isempty(colDepth) == 0
    sonic.Depth = u2u(1, inpDepthUnits, sonic.DepthUnits, [])*num(:,colDepth);
else
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Empty column %g containing the depths \n', colDepth);
    fprintf('>>> PAUSE -- Continue? \n');   pause;
end;

if isempty(colVp) == 0   
    tmp = num(:,colVp);
    if flagNaN == 1;   
        tmp(abs(tmp - missingData) < tol) = NaN;                % replace missing data with NaNs
    end;
    sonic.Vp  = u2u(1, inpVelocityLengthUnits, sonic.VelLengthUnits, [])/ ...
                u2u(1, inpVelocityTimeUnits,   sonic.VelTimeUnits,   [])*tmp;
elseif isempty(colQp) == 0
    tmp = num(:,colQp);
    if flagNaN == 1;   
        tmp(abs(tmp - missingData) < tol) = NaN;                % replace missing data with NaNs
    end;
    sonic.Vp  = u2u(1, inpSlownessLengthUnits, sonic.VelLengthUnits, [])/ ...
                u2u(1, inpSlownessTimeUnits,   sonic.VelTimeUnits,   [])./tmp;
else
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Empty columns containing the P-wave velocities and slownesses \n');
    fprintf('>>> PAUSE -- Continue? \n');   pause;
end;

if isempty(colVs1) == 0
    tmp = num(:,colVs1);
    if flagNaN == 1;   
        tmp(abs(tmp - missingData) < tol) = NaN;                % replace missing data with NaNs
    end;
    sonic.Vs1 = u2u(1, inpVelocityLengthUnits, sonic.VelLengthUnits, [])/ ...
                u2u(1, inpVelocityTimeUnits,   sonic.VelTimeUnits,   [])*tmp;
elseif isempty(colQs1) == 0
    tmp = num(:,colQs1);
    if flagNaN == 1;   
        tmp(abs(tmp - missingData) < tol) = NaN;                % replace missing data with NaNs
    end;
    sonic.Vs1 = u2u(1, inpSlownessLengthUnits, sonic.VelLengthUnits, [])/ ...
                u2u(1, inpSlownessTimeUnits,   sonic.VelTimeUnits,   [])./tmp;
else
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Empty columns containing the S1-wave velocities and slownesses \n');
    fprintf('>>> PAUSE -- Continue? \n');   pause;
end;

if isempty(colVs2) == 0
    tmp = num(:,colVs2);
    if flagNaN == 1;   
        tmp(abs(tmp - missingData) < tol) = NaN;                % replace missing data with NaNs
    end;
    sonic.Vs2 = u2u(1, inpVelocityLengthUnits, sonic.VelLengthUnits, [])/ ...
                u2u(1, inpVelocityTimeUnits,   sonic.VelTimeUnits,   [])*tmp;
elseif isempty(colQs2) == 0
    tmp = num(:,colQs2);
    if flagNaN == 1;   
        tmp(abs(tmp - missingData) < tol) = NaN;                % replace missing data with NaNs
    end;
    sonic.Vs2 = u2u(1, inpSlownessLengthUnits, sonic.VelLengthUnits, [])/ ...
                u2u(1, inpSlownessTimeUnits,   sonic.VelTimeUnits,   [])./tmp;
else
    % do nothing
end;

end    % of the function 