%% Tool for processing text files with ultrasonic experiment data
%  Geser Dugarov, 2018
%% Set paths and directories
clear variables; %close all;
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
fprintf('    is running... \n');
addpath(genpath('FuncBase'));

%% Settings
% path to the folder with data
path = {'c:\Geser\Dropbox\Projects\����������� - ��������������\���������� ��������\2018.11.15 - ��������� � ���\20181115\'};
% cell array of IDs
IDs = {'Al','Al78'};
paths = IDs;
for i = 1:numel(IDs)
    paths{i} = [path{1},IDs{i},'\'];
end
numLines = 10053; % number of lines in text files
skipLines = 3;   % number of lines that should be skipped
step = 1;        % step for skipping files (1 - without skip)
sett = struct(...
    'intWind',  0.5e-6, ... % s, window length for interpolatation and finding max
    'polyfitdegr', 3, ...      % degree of interpolation polynom
    'polyfitdisc', 1.0e-9, ... % interpolation discretization
    'noiseFilter', true,                    ... % filter noise or not
    'filtFreqRange', [0.0 0.0 1.5e6 1.6e6]  ... % filter before processing ./^^\.
    );

valuenum = numLines-skipLines;
meanData = nan(valuenum,numel(paths));
meanTime = nan(valuenum,numel(paths));
for pathID = 1:numel(paths)
    path = paths{pathID};
    ID = IDs{pathID};
    F = dir([path,'*.txt']); %����� ������ �� ��������������
    fileID = 1:step:length(F);
    tracenum = numel(fileID);
    Time  = nan(valuenum,tracenum);
    DataP = nan(valuenum,tracenum);
    DataS = nan(valuenum,tracenum);
    startPos = zeros(1,tracenum);
    count = 1;
    for i = fileID
        fname = F(i).name;
        temp = importdata([path,fname],'\t',skipLines);
        temp = temp.data;
        Time(:,count)  = temp(:,1);   
        DataP(:,count) = temp(:,3);
        DataS(:,count) = temp(:,6); 
        startPos(count) = find(Time(:,count) > 0,1);
        count = count+1;
    end

    if max(startPos)>mean(startPos) || min(startPos)<mean(startPos)
        error('Different length of traces');
    end

%     % remove data before 0 s
%     Time = Time(mean(startPos):end,:);
%     Data = Data(mean(startPos):end,:);
%     valuenum = valuenum-mean(startPos)+1;

    %% Change NaNs to 0
    DataP(isnan(DataP)) = 0;
    DataS(isnan(DataS)) = 0;
    %% Correcting shift up
    DataP = DataP - repmat(mean(DataP),size(DataP,1),1);
    DataS = DataS - repmat(mean(DataS),size(DataS,1),1);
    %% Remove some noise near 0 s
	DataP(mean(startPos):mean(startPos)+100,:) = 0;
	DataS(mean(startPos):mean(startPos)+100,:) = 0;
    
    %% Noise filtering
    if sett.noiseFilter
        for trace = 1:tracenum
            df = 1/(Time(end,trace)-Time(1,trace));
            freq = 0:df:df*(valuenum-1);
            filt = ones(1,size(freq,2));
            pos1 = floor(sett.filtFreqRange(1,1)/df)+1;
            pos2 = ceil(sett.filtFreqRange(1,2)/df)+1;
            pos3 = floor(sett.filtFreqRange(1,3)/df)+1;
            pos4 = ceil(sett.filtFreqRange(1,4)/df)+1;
            filt(1,[1:pos1 pos4:numel(filt)]) = 0;
            if (pos2-pos1-1 > 0)
                filt(1,pos1+1:pos2-1) = 1/(pos2-pos1):1/(pos2-pos1):(pos2-pos1-1)/(pos2-pos1);
            end
            if (pos4-pos3-1 > 0)
                filt(1,pos3+1:pos4-1) = (pos4-pos3-1)/(pos4-pos3):-1/(pos4-pos3):1/(pos4-pos3);
            end
            temp = fft(DataP(:,trace));
            temp = real(ifft(temp.*filt'));
            DataP(:,trace) = temp*(max(DataP(:,trace))/max(temp));
            temp = fft(DataS(:,trace));
            temp = real(ifft(temp.*filt'));
            DataS(:,trace) = temp*(max(DataS(:,trace))/max(temp));
            clear df freq filt pos1 pos2 pos3 pos4 temp;
        end
    end

    meanTime(:,pathID)  = mean(Time,2);
    meanDataP(:,pathID) = mean(DataP,2);
    meanDataS(:,pathID) = mean(DataS,2);
    
%     %% Plotting traces
%     FigTrace = figure();
%     set(FigTrace,'units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
%     subplot('Position',[0.05 0.06 0.90 0.88]); hold on;
%     % max values for normalization of traces
%     maxval = repmat(max(max(abs(Data))),1,size(Data,2));
%     tempData = repmat(1:size(Data,2),size(Data,1),1) + Data./repmat(maxval,size(Data,1),1)/2;
%     plot(Time,tempData,'k');
%     % tempData = (1:size(maxpos,1))' + maxpos(:,2)./maxval'/2;
%     % plot(maxpos(:,1),tempData,'+r');
%     % pos = find(maxpos(:,3));
%     % plot(maxpos(pos,1),tempData(pos,1),'+g');
%     clear maxval tempData pos;
%     xlim([0,Time(end,1)]);
%     ylim([0,20]);
%     pickspath = [path,'procTextData',ID,'.txt'];
%     LinePlotExplorer_polyfit(pickspath);
%     % figure settings for traces
%     set(gca,'Ytick',1:tracenum); grid on; xlabel('times, s'); ylabel('trace number'); title(ID,'Interpreter','none');
%     clear pickspath;

end

for wavetype = 1:2
    %% Load manual picks
    if wavetype < 2
        pickspath = 'data\fbs_hydrate_exp\manualPick_procTextData_meanDataP.dat';
        Data = DataP;
        meanData = meanDataP;
        waveID = 'P';
    else
        pickspath = 'data\fbs_hydrate_exp\manualPick_procTextData_meanDataS.dat';
        Data = DataS;
        meanData = meanDataS;
        waveID = 'S';
    end
    fileID = fopen(pickspath,'r');
    if fileID < 0
        manualPick = [];
    else
        manualPick = fscanf(fileID,'%d %f',[2 Inf]);
        manualPick = manualPick';
        fclose(fileID);
    end

    %% Find phase positions using the windows
    maxpos = nan(tracenum,3);
    if ~isempty(manualPick)
        for trace = 1:tracenum
            % accepting manual picks
            pos = find(manualPick(:,1) == trace);
            if ~isempty(pos)
                temptime = manualPick(pos,2);
                % cut signal using intWind and interpolate
                cutSig = cutSignal([meanTime(:,trace) meanData(:,trace)], ...
                                   [temptime-sett.intWind/2,temptime+sett.intWind/2]);
                inttime = cutSig(1,1):sett.polyfitdisc:cutSig(end,1);
                ws = warning('off','all');
                intpoly = polyfit(cutSig(:,1),cutSig(:,2),sett.polyfitdegr);
                warning(ws);
                [maxpos(trace,3),pos] = max(abs(polyval(intpoly,inttime)));
                if polyval(intpoly,inttime(1,pos))<0
                    maxpos(trace,3) = -maxpos(trace,3);
                end
                maxpos(trace,2) = inttime(1,pos);
                maxpos(trace,1) = trace;
            end
        end
    end
    maxpos(isnan(maxpos(:,1)),:) = [];
    clear pos temptime inttime intpoly trace ws cutSig;

    %% Plotting traces
    FigTrace = figure();
    set(FigTrace,'units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
    subplot('Position',[0.05 0.06 0.90 0.88]); hold on;
    % max values for normalization of traces
    normFlag = true;
    if normFlag
        maxval = max(abs(meanData));
        normStr = 'normalized';
    else
        maxval = repmat(max(max(abs(meanData))),1,size(meanData,2));
        normStr = 'true';
    end
    tempData = repmat(1:size(meanData,2),size(meanData,1),1) + meanData./repmat(maxval,size(meanData,1),1)/2;
    plot(meanTime,tempData,'k');
    if ~isempty(maxpos)
        tempData = maxpos(:,1) + maxpos(:,3)./maxval(maxpos(:,1))'/2;
        plot(maxpos(:,2),tempData,'+r');
    end
    clear maxval tempData pos;
    xlim([0,meanTime(end,1)]);
    ylim([0,6.9]);
    LinePlotExplorer_polyfit(pickspath);
    % figure settings for traces
    set(gca,'Ytick',1:tracenum); grid on; xlabel('times, s'); ylabel('trace number');
    title(['Mean data, ',normStr,' amplitudes, ',waveID,' wave']);
    clear pickspath;

end

fprintf('>>> Return from script ''%s'' \n \n', [thisFileName, '.m']);
clear thisFileName thisFolderName;
return