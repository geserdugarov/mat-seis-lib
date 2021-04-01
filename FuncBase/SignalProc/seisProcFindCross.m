%% *seisProcFindCross*
% Find zero crosses

%% 
% *Input:* 
% Time      - [N]            time array
% Data      - [N, traceNum]  traces
% maxpos    - [traceNum, 3]  first break data

%%
% *Output:*
% cross0 - [traceNum, 2] zero cross data

%%
% *Author:* Geser Dugarov 2018

%%

function [cross0] = seisProcFindCross(Time, Data, maxpos)

% [thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));

tracenum = size(Data,2);
cross0  = NaN(tracenum,2);
for trace = 1:tracenum
    data = Data(:,trace);
    time = Time;
    data = data(time<maxpos(trace,1));
    time = time(time<maxpos(trace,1));
    [~,t0] = crossing(data,time,0,'linear');
    cross0(trace,1) = t0(end-1);
    data = Data(:,trace);
    time = Time;
    data = data(time>maxpos(trace,1));
    time = time(time>maxpos(trace,1));
    [~,t0] = crossing(data,time,0,'linear');
    cross0(trace,2) = t0(2);
end

end