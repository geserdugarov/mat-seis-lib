%% *seisProcCalcStd*
% The amplitude std estimation

%% 
% *Input:* 
% Time      - [N]            time array
% Data      - [N, traceNum]  traces
% cross0    - [traceNum, 1]  zero cross data

%%
% *Output:*
% resRMS    - [traceNum] estimated std of amplitudes

%%
% *Author:* Geser Dugarov 2018

%%

function [ener] = seisProcCalcStd(Time, Data, cross0)

% [thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));

tracenum = size(Data,2);
ener = NaN(tracenum,1);
for trace = 1:tracenum
    data = Data(:,trace);
    data = data( Time>cross0(trace,1) & Time<cross0(trace,2) );
    ener(trace,1) = std(data);
end

end