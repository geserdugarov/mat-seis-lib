%% *windows*
% Adding window to signal

%% 
% *Input:* 
% trace - [times (s), amplitudes] seismogram trace
% dtw   - [scalar (s)] length of smoothing
% type  - [string] type of smoothing: barthannwin, bartlett,
%                  blackman or hamming

%%
% *Output:*
% sigCut - [times (s), amplitudes] cutted windowed signal

%%
% *Author:* Geser Dugarov 2016

%%
function [trace] = windows(trace, dtw, type)

[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));

if (dtw > trace(end,1)-trace(1,1))
    error('>>> Smoothing window > signal length.\n');
end

if (dtw>0)
    timestep = trace(2,1) - trace(1,1);
    pointNumForWind = round(dtw/timestep);
    switch type
        case 'barthannwin'
            window = barthannwin(2*pointNumForWind);
        case 'bartlett'
            window = bartlett(2*pointNumForWind);
        case 'blackman'
            window = blackman(2*pointNumForWind);
        case 'hamming'
            window = hamming(2*pointNumForWind);
        otherwise
            fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
            fprintf('>>> There is no window type: ''%s''\n\n',type);
            error('>>> STOP');
    end
	trace(1:pointNumForWind+1, 2) = trace(1:pointNumForWind+1, 2).*[0; window(1:pointNumForWind)];
	trace(size(trace,1)-pointNumForWind:size(trace,1), 2) = trace(size(trace,1)-pointNumForWind:size(trace,1), 2).*[window(pointNumForWind+1:2*pointNumForWind); 0];
end

end % of the function
