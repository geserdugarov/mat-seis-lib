%% *seisProcCalcEner*
% The wave energy estimation

%% 
% *Input:* 
% Time      - [N]            time array
% Data      - [N, traceNum]  traces
% cross0    - [traceNum, 1]  zero cross data

%%
% *Output:*
% ener      - [traceNum] estimated energy

%%
% *Author:* Geser Dugarov 2018

%%

function [ener] = seisProcCalcEner(Time, Data, cross0)

% [thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));

tracenum = size(Data,2);
ener = NaN(tracenum,1);
for trace = 1:tracenum
    data = Data(:,trace);
    data = data( Time>cross0(trace,1) & Time<cross0(trace,2) );
    ener(trace,1) = sum(data.^2);
end

end