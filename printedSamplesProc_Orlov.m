%% Tool for processing ultrasonic experiment data
%  Velocity estimation
%  Max-phase signal is expected
%  Geser Dugarov, 2019

%% Set paths and directories
clear variables; close all;
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
fprintf('    is running... \n');
addpath(genpath('FuncBase'));

%% Data director path
datadir = 'd:\Dropbox\Projects\грант 2019-2021 РФФИ а (3D-печать трещиноватости)\акустическое просвечивание\';

filepath = {'P_M1_M5\1','M1\4','M1\1'};
fileIDs  = {'P','S_1','S_2'};
ID = 'FDM, изотропный'; pickpath = [datadir,'manualpick_',ID];
smplsize = [36;36;36];
f2f = [1.2e-6; 2.25e-6; 2.25e-6];

Time = repmat((0:0.5e-8:9.9995e-5)',1,numel(filepath));
gain = 1.5;
tracenum = numel(filepath);

%% Settings for V estimation
settV = struct(...
    'useWind',     true, ...
    'searchWind',  [1.5e-5 2.0e-5; 3.4e-5 4.0e-5; 3.6e-5 4.0e-5], ...
    'firstVal', 0.45, ...    % search first settV.firstVal*max 0.45
    'maxWind',  0.2e-5, ... % s, window after range*max
    'intWind',  0.06e-5, ... % s, window length for interpolatation and finding max
    'polyfitdegr', 3, ...      % degree of interpolation polynom
    'polyfitdisc', 1.0e-9, ... % interpolation discretization
    'noiseFilter', true,   ... % filter noise or not
    'filtFreqRange', [0 50e3 1e6 2e6;  0 50e3 350e3 500e3;  0 50e3 350e3 500e3], ... % ,,/''\,,
    'useF2Fmean', false, ... % use mean correction or curent f2f data for sample
    'plotTracesNorm', true ... % true - plot normalized traces, false - real amplitudes
    );

%% Processing
type = 'PS';
Velocities = zeros(tracenum,1);
    
for trace=1:tracenum
    Data(:,trace) = importdata([datadir,filepath{trace}],'\t',0);
end
% Change NaNs to 0
Data(isnan(Data)) = 0;
% Correcting shift up
Data = Data - repmat(mean(Data),size(Data,1),1);

% Load manual picks for current identif
fileID = fopen(pickpath,'r');
if fileID < 0
    manualPick = [];
else
    manualPick = fscanf(fileID,'%d %f',[2 Inf]);
    manualPick = manualPick';
    fclose(fileID);
end

%% Noise filtering
if settV.noiseFilter
    for trace = 1:tracenum
        % choose time (f2f, aluminum or sample)
        df = 1/(Time(end,trace)-Time(1,trace));
        freq = 0:df:df*(size(Time,1)-1);
        filt = ones(1,size(freq,2));
        pos1 = floor(settV.filtFreqRange(trace,1)/df)+1;
        pos2 = ceil(settV.filtFreqRange(trace,2)/df)+1;
        pos3 = floor(settV.filtFreqRange(trace,3)/df)+1;
        pos4 = ceil(settV.filtFreqRange(trace,4)/df)+1;
        filt(1,[1:pos1 pos4:numel(filt)]) = 0;
        if (pos2-pos1-1 > 0)
            filt(1,pos1+1:pos2-1) = 1/(pos2-pos1):1/(pos2-pos1):(pos2-pos1-1)/(pos2-pos1);
        end
        if (pos4-pos3-1 > 0)
            filt(1,pos3+1:pos4-1) = (pos4-pos3-1)/(pos4-pos3):-1/(pos4-pos3):1/(pos4-pos3);
        end
        temp = fft(Data(:,trace));
        temp = real(ifft(temp.*filt'));
        Data(:,trace) = temp*(max(Data(:,trace))/max(temp));
        clear df freq filt pos1 pos2 pos3 pos4 temp;
    end
end

%% Find first max position using the windows
for trace = 1:tracenum
    % accepting manual picks
    temptime = NaN;
    if ~isnan(manualPick)
        pos = find(manualPick(:,1) == trace);
        if ~isempty(pos)
            temptime = manualPick(pos,2);
            maxpos(trace,3) = 1;
        end
    end
    % if there is no manual pick
    if isnan(temptime)
        if settV.useWind % search in a window
            time = Time( Time(:,trace)>settV.searchWind(trace,1) & Time(:,trace)<settV.searchWind(trace,2),trace);
            data = Data( Time(:,trace)>settV.searchWind(trace,1) & Time(:,trace)<settV.searchWind(trace,2),trace);
            posWind = sum(Time(:,trace) <= settV.searchWind(trace,1));
        else % search in a whole time range
            time = Time(:,trace);
            data = Data(:,trace);
            posWind = 0;
        end
        % find max positions for polynomial interpolation
        findval = settV.firstVal*max(data);
        rangepos = find(data>findval,1,'first');
        if isempty(rangepos)
            maxpos(trace,:) = [NaN NaN NaN];
            continue;
        end
        % cut signal using maxWind and find max
        cutSig = cutSignal([time data], ...
                           [time(rangepos),time(rangepos)+settV.maxWind]);
        [val,pos] = max(cutSig(:,2));
        temptime = Time(posWind+rangepos+pos-1,trace);
        maxpos(trace,3) = 0;
    end
    % cut signal using intWind and interpolate
    cutSig = cutSignal([Time(:,trace) Data(:,trace)], ...
                       [temptime-settV.intWind/2,temptime+settV.intWind/2]);
    inttime = cutSig(1,1):settV.polyfitdisc:cutSig(end,1);
    ws = warning('off','all');
    intpoly = polyfit(cutSig(:,1),cutSig(:,2),settV.polyfitdegr);
    warning(ws);
    [maxpos(trace,2),pos] = max(polyval(intpoly,inttime));
    maxpos(trace,1) = inttime(1,pos);

    %% Velocity calculation
    Velocities(trace) = 1.0e-3*smplsize(trace)/(maxpos(trace,1)-f2f(trace,1)); % mm/s -> m/s
end

FigTrace = figure;
set(FigTrace,'units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
subplot('Position',[0.03 0.10 0.95 0.85]); hold on;
if settV.plotTracesNorm
    maxval = max(abs(Data));
else
    maxval = repmat(max(max(abs(Data))),1,tracenum);
end
tempData = repmat(1:size(Data,2),size(Data,1),1) + Data./repmat(maxval,size(Data,1),1)/2*gain;
plot(Time, tempData, 'k');

        % markers for first brakes (by phase)
        tempData = (1:size(maxpos,1))' + maxpos(:,2)./maxval'/2*gain;
        plot(maxpos(:,1),tempData,'+r');
        pos = find(maxpos(:,3));
        plot(maxpos(pos,1),tempData(pos,1),'+g');

xlim([0,max(Time(end,:))]); ylim([0,tracenum+1]);
LinePlotExplorer_polyfit(pickpath);
set(gca,'Ytick',1:tracenum);
grid on;
xlabel('Время, с');
title(ID);
set(gca,'YTickMode','manual'); set(gca,'YTick',1:tracenum);
set(gca,'YTickLabelMode','manual'); set(gca,'YTickLabel',fileIDs);
set(gca,'FontSize',16);


fprintf('>>> Return from script ''%s'' \n \n', [thisFileName, '.m']);
clear thisFileName thisFolderName;
return