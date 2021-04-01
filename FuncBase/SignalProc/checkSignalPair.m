%% *checkSignalPair*
% Check length and discretization of two signals

%% 
% *Input:* 
% trace1 - [times, amplitudes] seismogram trace
% trace2 - [times, amplitudes] seismogram trace

%%
% *Output:*
% trace1 - [times, amplitudes] seismogram trace
% trace2 - [times, amplitudes] seismogram trace

%%
% *Author:* Geser Dugarov 2016

%%
function [trace1, trace2] = checkSignalPair(trace1, trace2)

[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));

checkSignalAndDiscr(trace1);
checkSignalAndDiscr(trace2);

t1 = trace1(end,1)-trace1(1,1);
t2 = trace2(end,1)-trace2(1,1);
if (100*abs(t1-t2)/mean([t1 t2]) > 0.1)% if time length differences more than 0.1% of mean
	fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> WARNING: Not equal time length of signals. \n');
    fprintf('>>> WARNING: Skipping some points. \n \n');
    dt1 = trace1(2,1)-trace1(1,1);
    dt2 = trace2(2,1)-trace2(1,1);
    if (rem(max([dt1 dt2]),min([dt1 dt2])) < 1.0e-6)
        temp = [];
        if (t1 > t2)
            for i=1:size(trace1,1)
                if (abs(trace1(i,1)-trace1(1,1)) <= t2)
                    temp = [temp; trace1(i,:)];
                end
            end
            trace1 = temp;
        else
            for i=1:size(trace2,1)
                if (abs(trace2(i,1)-trace2(1,1)) <= t1)
                    temp = [temp; trace2(i,:)];
                end
            end
            trace2 = temp;
        end
    else
        fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
        fprintf('>>> Different discretization. Can not skip some points. \n \n');
        error('>>> STOP');
    end
end

end % of the function
