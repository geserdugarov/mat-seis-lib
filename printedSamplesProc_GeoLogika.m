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
datadir = 'd:\Dropbox\Projects\грант 2019-2021 РФФИ а (3D-печать трещиноватости)\3D model\acoustic\2019.10.18 - GeoLogika\';

% filenames = {'f2f','sH00', 'sH45', 'sH90', 'sV00', ...
%                    'e1H00','e1H45','e1H90','e1V00', ...
%                    'e2H00','e2H45','e2H90','e2V00'};
% pickpath = [datadir,'manualpick_1_'];
% smplsize = repmat([repmat(36.2,3,1); 40],3,1);

filenames = {'f2f','sH00', 'sH45', 'sH90', 'sV00', 'sV45', 'sV90', ...
                   'e1H00','e1H45','e1H90','e1V00','e1V45','e1V90', ...
                   'e2H00','e2H45','e2H90','e2V00','e2V45','e2V90'};
pickpath = [datadir,'manualpick_2_'];
smplsize = repmat([repmat(36.2,3,1); repmat(40,3,1)],3,1);

tracenum = numel(filenames);

%% Settings for V estimation
settV = struct(...
    'useWind',     true, ...
    'searchWind',  [30 50; 60 90], ...
    'firstVal', 0.45, ...    % search first settV.firstVal*max 0.45
    'maxWind',  3, ... % mks, window after range*max
    'intWind',  1.5, ... % mks, window length for interpolatation and finding max
    'polyfitdegr', 3, ...      % degree of interpolation polynom
    'polyfitdisc', 1.0e-2, ... % interpolation discretization
    'noiseFilter', true,                    ... % filter noise or not
    'filtFreqRange', [0.0 0.0 1.5e6 1.6e6], ... % filter before processing ./^^\.
    'useF2Fmean', false, ... % use mean correction or curent f2f data for sample
    'plotTracesNorm', true ... % true - plot normalized traces, false - real amplitudes
    );

%% Processing
type = 'PS';
Velocities = zeros(tracenum-1,2);
for wavetype = 1:2
    for trace=1:tracenum
        data = dlmread([datadir,filenames{trace},'.csv'],';',34,0);
        Time(:,trace) = data(:,1);
        Data(:,trace) = data(:,wavetype*3-1);
    end
    % Change NaNs to 0
    Data(isnan(Data)) = 0;
    % Correcting shift up
    Data = Data - repmat(mean(Data),size(Data,1),1);
    
    % Load manual picks for current identif
    fileID = fopen([pickpath,type(wavetype),'.txt'],'r');
    if fileID < 0
        manualPick = [];
    else
        manualPick = fscanf(fileID,'%d %f',[2 Inf]);
        manualPick = manualPick';
        fclose(fileID);
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
                    time = Time( Time(:,trace)>settV.searchWind(wavetype,1) & Time(:,trace)<settV.searchWind(wavetype,2),trace);
                    data = Data( Time(:,trace)>settV.searchWind(wavetype,1) & Time(:,trace)<settV.searchWind(wavetype,2),trace);
                    posWind = sum(Time(:,trace) <= settV.searchWind(wavetype,1));
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
        end

        %% Velocity calculation
        Velocities(:,wavetype) = 1.0e-3*smplsize/1.0e-6./(maxpos(2:end,1)-maxpos(1,1)); % mm/mks -> m/s
        
    FigTrace = figure;
    set(FigTrace,'units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
    subplot('Position',[0.05 0.06 0.90 0.88]); hold on;
    if settV.plotTracesNorm
        maxval = max(abs(Data));
    else
        maxval = [max(abs(Data(:,1))) max(abs(Data(:,2))) ...
                  repmat(max(max(abs(Data(:,3:end)))),1,size(Data,2)-2)];
    end
    tempData = repmat(1:size(Data,2),size(Data,1),1) + Data./repmat(maxval,size(Data,1),1)/2;
    plot(Time, tempData, 'k');
    
            tempData = (1:size(maxpos,1))' + maxpos(:,2)./maxval'/2;
            plot(maxpos(:,1),tempData,'+r');
            pos = find(maxpos(:,3));
            plot(maxpos(pos,1),tempData(pos,1),'+g');
            
    xlim([0,max(Time(end,:))]); ylim([0,tracenum+1]);
    LinePlotExplorer_polyfit([pickpath,type(wavetype),'.txt']);
    set(gca,'Ytick',1:tracenum);
    grid on;
    xlabel('times, \mus');
    title([type(wavetype),' waves']);
    set(gca,'YTickMode','manual'); set(gca,'YTick',1:tracenum);
    set(gca,'YTickLabelMode','manual'); set(gca,'YTickLabel',filenames);
end



fprintf('>>> Return from script ''%s'' \n \n', [thisFileName, '.m']);
clear thisFileName thisFolderName;
return