%% *cutSignal*
% Signal cutting with specified window

%% 
% *Input:* 
% trace - [times (s), amplitudes] seismogram trace
% twind    - [scalar (s) scalar (s)] time window

%%
% *Output:*
% sigCut - [times (s), amplitudes] cutted windowed signal

%%
% *Author:* Geser Dugarov 2016

%%
function [sigCut] = cutSignal(trace, twind)

[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));

time = trace(:,1);
data = trace(:,2);
data = data(time>=twind(1));
time = time(time>=twind(1));
data = data(time<=twind(2));
time = time(time<=twind(2));
sigCut = [time data];

if (isempty(sigCut))
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
	error('>>> There is no signal after cutting, please check input parameters.');
end

end % of the function
