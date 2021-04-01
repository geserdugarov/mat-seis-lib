%% *bandpassZeroPhase*
% Bandpass zero-phase filtering

%%
% *Input:*

% x         - [:, 1] time series
% dTime     - [scalar] time sample  
% freqBand  - [1, 4] array of four frequencies defining a trapezoidal filter
% timeUnits - [string] units of time pertaining to variable dTime
% freqUnits - [string] units of frequency pertaining to array freqBand

%%
% *Output:*

% y         - filtered time series

%%
% *Author:* Vladimir Grechka 2012

%%
function [y] = bandpassZeroPhase(x, dTime, freqBand, timeUnits, freqUnits)
%% Settings and defaults
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

dTimeS = dTime;   freqBandHz = freqBand;

if nargin < 5  ||  isempty(freqUnits) == 1
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Unspecified frequency units defaulted to [Hz] -- PAUSE \n \n');   pause;
else
    freqBandHz = u2u(freqBand, freqUnits, 'Hz', 1); 
end;

if nargin < 4  ||  isempty(timeUnits) == 1
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Unspecified time units defaulted to [s] -- PAUSE \n \n');   pause;
else
    dTimeS = u2u(dTime, timeUnits, 's', 1); 
end;

if length(freqBand) ~= 4
    fprintf('\n>>> Function ''%s'' \n', functionInfo.file);
    fprintf('>>> Length of trapezoidal filter ''freqBand'' should be 4 rather than %g \n \n ', ...
            length(freqBand));
      error('>>> STOP');
end;

%% Design recursive digital filter 
nyquistHz = 0.5/dTimeS;       % Nyquist frequency in Hz
[b, a] = yulewalk(10, [0, freqBandHz/nyquistHz, 1], [0, 0, 1, 1, 0, 0]);

% Zero-phase filtering
y = filtfilt(b, a, detrend(x));   

end    % of the function