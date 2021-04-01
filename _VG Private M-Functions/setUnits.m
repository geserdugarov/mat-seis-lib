%% *setUnits*
% Set up structure containing units of length, time, and angle to be used for data input, processing, 
% and displaying the results

%%
% *Input:* 

% unitsLengthInput   - [string] input units of length
% unitsLengthProc    - [string] units of length used in data processing
% unitsLengthDisplay - [string] units of length used for displays
% unitsTimeInput     - [string] input units of time
% unitsTimeProc      - [string] units of time used in data processing
% unitsTimeDisplay   - [string] units of time used for displays
% unitsAngleInput    - [string] input units of angle
% unitsAngleProc     - [string] units of angle used in data processing
% unitsAngleDisplay  - [string] units of angle used for displays
% unitsFreqInput     - [string] input units of frequency
% unitsFreqProc      - [string] units of frequency used in data processing
% unitsFreqDisplay   - [string] units of frequency used for displays

%%
% *Output:* 

% units.             - nested structure containing the following fields:
%     length.
%         input      = 'unitsLengthInput'
%         processing = 'unitsLengthProc'
%         display    = 'unitsLengthDisplay'
%     time.
%         input      = 'unitsTimeInput'
%         processing = 'unitsTimeProc'
%         display    = 'unitsTimeDisplay'
%     angle.
%         input      = 'unitsAngleInput'
%         processing = 'unitsAngleProc'
%         display    = 'unitsAngleDisplay'
%     frequency.
%         input      = 'unitsFreqInput'
%         processing = 'unitsFreqProc'
%         display    = 'unitsFreqDisplay'

%% 
% *Author:* Vladimir Grechka 2012

function [units] = setUnits(unitsLengthInput, unitsLengthProc, unitsLengthDisplay, ...
                            unitsTimeInput,   unitsTimeProc,   unitsTimeDisplay,   ...
                            unitsAngleInput,  unitsAngleProc,  unitsAngleDisplay,  ...
                            unitsFreqInput,   unitsFreqProc,   unitsFreqDisplay)
%% Settings and defaults
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

if isempty(unitsLengthInput) == 1
    unitsLengthInput = 'm';                       % default for variable 'unitsLengthInput'
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Empty input variable ''unitsLengthInput'' defaulted to ''%s'' \n', unitsLengthInput);
end;

if isempty(unitsLengthProc) == 1
    unitsLengthProc = unitsLengthInput;           % default for variable 'unitsLengthProc'
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Empty input variable ''unitsLengthProc'' defaulted to ''%s'' \n', unitsLengthProc);
end;

if isempty(unitsLengthDisplay) == 1
    unitsLengthDisplay = unitsLengthInput;        % default for variable 'unitsLengthDisplay'
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Empty input variable ''unitsLengthDisplay'' defaulted to ''%s'' \n', unitsLengthDisplay);
end;

if isempty(unitsTimeInput) == 1
    unitsTimeInput = 's';                         % default for variable 'unitsTimeInput'
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Empty input variable ''unitsTimeInput'' defaulted to ''%s'' \n', unitsTimeInput);
end;

if isempty(unitsTimeProc) == 1
    unitsTimeProc = unitsTimeInput;               % default for variable 'unitsTimeProc'
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Empty input variable ''unitsTimeProc'' defaulted to ''%s'' \n', unitsTimeProc);
end;

if isempty(unitsTimeDisplay) == 1
    unitsTimeDisplay = unitsTimeInput;            % default for variable 'unitsTimeDisplay'
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Empty input variable ''unitsTimeDisplay'' defaulted to ''%s'' \n', unitsTimeDisplay);
end;

if isempty(unitsAngleInput) == 1
    unitsAngleInput = 'rad';                      % default for variable 'unitsAngleInput'
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Empty input variable ''unitsAngleInput'' defaulted to ''%s'' \n', unitsAngleInput);
end;

if isempty(unitsAngleProc) == 1
    unitsAngleProc = unitsAngleInput;             % default for variable 'unitsAngleProc'
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Empty input variable ''unitsAngleProc'' defaulted to ''%s'' \n', unitsAngleProc);
end;

if isempty(unitsAngleDisplay) == 1
    unitsAngleDisplay = unitsAngleInput;          % default for variable 'unitsAngleDisplay'
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Empty input variable ''unitsAngleDisplay'' defaulted to ''%s'' \n', unitsAngleDisplay);
end;

if isempty(unitsFreqInput) == 1
    unitsFreqInput = 'Hz';                        % default for variable 'unitsFreqInput'
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Empty input variable ''unitsFreqInput'' defaulted to ''%s'' \n', unitsFreqInput);
end;

if isempty(unitsFreqProc) == 1
    unitsFreqProc = unitsFreqInput;               % default for variable 'unitsFreqProc'
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Empty input variable ''unitsFreqProc'' defaulted to ''%s'' \n', unitsFreqProc);
end;

if isempty(unitsFreqDisplay) == 1
    unitsFreqDisplay = unitsFreqInput;            % default for variable 'unitsFreqDisplay'
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Empty input variable ''unitsFreqDisplay'' defaulted to ''%s'' \n', unitsFreqDisplay);
end;


%% Fill out the output 'units' structure
units.length.input         = unitsLengthInput;
units.length.processing    = unitsLengthProc;
units.length.display       = unitsLengthDisplay;
units.time.input           = unitsTimeInput;
units.time.processing      = unitsTimeProc;
units.time.display         = unitsTimeDisplay;
units.angle.input          = unitsAngleInput;
units.angle.processing     = unitsAngleProc;
units.angle.display        = unitsAngleDisplay;
units.frequency.input      = unitsFreqInput;
units.frequency.processing = unitsFreqProc;
units.frequency.display    = unitsFreqDisplay;

end    % of the function