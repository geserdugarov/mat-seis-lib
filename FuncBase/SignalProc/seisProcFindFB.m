%% *seisProcFindFB*
% Find first breaks

%% 
% *Input:* 
% Time       - [N]              time array
% Data       - [N, traceNum]    traces
% sett       - [struct]         setting structure
% pickspath  - [string]         path to manual picks
% predefPick - [traceNum times] predefined times for max values

% sett = struct(...
%     'searchWind',  [0.5 0.8], ... % time windows for wave searching
%     'firstVal',    0.6, ...     % search first sett.firstVal*max
%     'firstValMod', true, ...    % if true, then max from |x|, else from x
%     'maxWind',  0.030, ... % s, window after range*max
%     'intWind',  0.007, ... % s, window length for interpolatation and finding max
%     'polyfitdegr', 3, ...      % degree of interpolation polynom
%     'polyfitdisc', 1.0e-5 ... % interpolation discretization
%     );

%%
% *Output:*
% maxpos - [Times (s), amplitudes, manual or auto] max position data

%%
% *Author:* Geser Dugarov 2018

%%

function [maxpos] = seisProcFindFB(Time, Data, sett, pickspath, predefPick)

% [thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));

% read data for manual picks from file
fileID = fopen(pickspath,'r');
if fileID < 0
    manualPick = [];
else
    manualPick = fscanf(fileID,'%d %f',[2 Inf]);
    manualPick = manualPick';
    fclose(fileID);
end
% processing predefined picks
if isempty(predefPick)
    manualPick = [];
else
    manualPick = predefPick;
end

% searching first break for each trace
tracenum = size(Data,2);
maxpos  = NaN(tracenum,3);
for trace = 1:tracenum
    % accepting manual picks
    % tempTime - Time for max|.|
    tempTime = NaN;
    if ~isnan(manualPick)
        pos = find(manualPick(:,1) == trace);
        if ~isempty(pos)
            tempTime = manualPick(pos,2);
            maxpos(trace,3) = 1;
        end
    end
    % if there is no manual pick
    if isnan(tempTime)
        % find max positions for polynomial interpolation
        if sett.useWind % search in a window
            time = Time( Time>sett.searchWind(1) & Time<sett.searchWind(2));
            data = Data( Time>sett.searchWind(1) & Time<sett.searchWind(2) ,trace);
            posWind = sum(Time <= sett.searchWind(1));
        else % search in a whole time range
            data = Data(:,trace);
            posWind = 0;
        end
        if sett.firstValMod % search using |x|
            findval = sett.firstVal*max(abs(data));
            rangepos = find(abs(data)>findval,1,'first');
        else % search using x
            findval = sett.firstVal*max(data);
            rangepos = find(data>findval,1,'first');
        end
        if isempty(rangepos)
            maxpos(trace,:) = [NaN NaN NaN];
            continue;
        end
        % cut signal using maxWind and find max
        cutSig = cutSignal([time data], ...
                           [time(rangepos),time(rangepos)+sett.maxWind]);
        if sett.firstValMod
            [~,pos] = max(abs(cutSig(:,2)));
        else
            [~,pos] = max(cutSig(:,2));
        end
        tempTime = Time(posWind+rangepos+pos-1);
        tempval  = Data(posWind+rangepos+pos-1,trace);
        maxpos(trace,3) = 0;
    else
        [~,pos] = min(abs(Time - tempTime));
        tempval = Data(pos,trace);
    end
    % cut signal using intWind and interpolate
    cutSig = cutSignal([Time Data(:,trace)], ...
                       [tempTime-sett.intWind/2,tempTime+sett.intWind/2]);
    intTime = cutSig(1,1):sett.polyfitdisc:cutSig(end,1);
    ws = warning('off','all');
    intpoly = polyfit(cutSig(:,1),cutSig(:,2),sett.polyfitdegr);
    warning(ws);
    if tempval > 0
        [maxpos(trace,2),pos] = max(polyval(intpoly,intTime));
    else
        [maxpos(trace,2),pos] = min(polyval(intpoly,intTime));
    end
    maxpos(trace,1) = intTime(1,pos);
end

end