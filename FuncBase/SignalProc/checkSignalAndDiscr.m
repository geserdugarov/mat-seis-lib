%% *checkSignalAndDiscr*
% Check signal and discretization

%% 
% *Input:* 
% trace - [times (s), amplitudes] seismogram trace

%%
% *Output:*
% no output

%%
% *Author:* Geser Dugarov 2016

%%
function [ ] = checkSignalAndDiscr(trace)

[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));

if (max(isnan(trace(:,1)))>0)
	fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> There are NaN time values. \n \n');
	error('>>> STOP');
end
if (max(isnan(trace(:,2)))>0)
	fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> There are NaN signal values. \n \n');
	error('>>> STOP');
end
dt = abs(trace(2:end,1)-trace(1:end-1,1));
if (100*(mean(dt)-max(dt))/mean(dt) > 0.1) % if deviation of mean from max more than 0.1%
	fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Bad discretization step, please check time values. \n \n');
	error('>>> STOP');
end

end % of the function
