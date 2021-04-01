clear all; close all;
addpath(genpath('FuncBase'));

sett = struct(...
    'intWind',  0.8e-6, ... % s, window length for interpolatation and finding max
    'polyfitdegr', 3, ...      % degree of interpolation polynom
    'polyfitdisc', 1.0e-9, ... % interpolation discretization
    'noiseFilter', true,                    ... % filter noise or not
    'filtFreqRange', [0.0 0.0 1.5e6 1.6e6]  ... % filter before processing ./^^\.
    );

datadir = 'c:\Geser\hydrates_data\';
filenamelist = dir([datadir,'*f2f*.mat']);
filenamelistID = false(numel(filenamelist),1);
for i = 1:size(filenamelist,1)
    name = filenamelist(i).name;
    if str2double(name(1:4))<2017
        filenamelistID(i) = true;
    else
        if str2double(name(1:4))==2017 && str2double(name(5:6))<7
            filenamelistID(i) = true;
        end
    end
end
filenamelist(filenamelistID) = [];

f2fTime  = [];
f2fDataP = [];
f2fDataS = [];
AlTime  = [];
AlDataP = [];
AlDataS = [];
for i = 1:numel(filenamelist)
    etaldata = load([datadir,filenamelist(i).name]);
    f2fTime  = [f2fTime  etaldata.Time(:,1)];
    f2fDataP = [f2fDataP etaldata.DataP(:,1)];
    f2fDataS = [f2fDataS etaldata.DataS(:,1)];
    AlTime  = [AlTime  etaldata.Time(:,1)];
    AlDataP = [AlDataP etaldata.DataP(:,7)];
    AlDataS = [AlDataS etaldata.DataS(:,7)];
end
clear etaldata i;
tracenum = numel(filenamelist);
valuenum = size(f2fTime,1);

%% Noise filtering
if sett.noiseFilter
    for i = 1:4
        switch i
            case 1
                Time = f2fTime;
                Data = f2fDataP;
            case 2
                Time = f2fTime;
                Data = f2fDataS;
            case 3
                Time = AlTime;
                Data = AlDataP;
            case 4
                Time = AlTime;
                Data = AlDataS;
        end
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
            temp = fft(Data(:,trace));
            temp = real(ifft(temp.*filt'));
            Data(:,trace) = temp*(max(Data(:,trace))/max(temp));
            clear df freq filt pos1 pos2 pos3 pos4 temp;
        end
    end
end

%% Find phase positions using the windows
maxpos = nan(tracenum,2,4);
for i = 1:4
    switch i
        case 1
            Time = f2fTime;
            Data = f2fDataP;
            temptime = 5.04e-6;
        case 2
            Time = f2fTime;
            Data = f2fDataS;
            temptime = 8.05e-6;
        case 3
            Time = AlTime;
            Data = AlDataP;
            temptime = 11.14e-6;
        case 4
            Time = AlTime;
            Data = AlDataS;
            temptime = 20.86e-6;
    end
    for trace = 1:tracenum
        % cut signal using intWind and interpolate
        cutSig = cutSignal([Time(:,trace) Data(:,trace)], ...
                           [temptime-sett.intWind/2,temptime+sett.intWind/2]);
        inttime = cutSig(1,1):sett.polyfitdisc:cutSig(end,1);
        ws = warning('off','all');
        intpoly = polyfit(cutSig(:,1),cutSig(:,2),sett.polyfitdegr);
        warning(ws);
        [maxpos(trace,2,i),pos] = max(abs(polyval(intpoly,inttime)));
        if polyval(intpoly,inttime(1,pos))<0
            maxpos(trace,2,i) = -maxpos(trace,2,i);
        end
        maxpos(trace,1,i) = inttime(1,pos);
    end
end
clear temptime inttime intpoly trace ws cutSig;

%% Estimation of velocities
Vp = 40.026/1000./(maxpos(:,1,3)-maxpos(:,1,1));
Vs = 40.026/1000./(maxpos(:,1,4)-maxpos(:,1,2));
fprintf('Mean Vp = %.0f\n',mean(Vp));
fprintf('Mean Vs = %.0f\n',mean(Vs));

%% Plot data
colorStep = 200/255/numel(filenamelist)*3;
temp1 = (1:-colorStep:0)';
zeroT = zeros(size(temp1));
colMap = [[temp1 zeroT zeroT]; ...
          [zeroT temp1 zeroT]; ...
          [zeroT zeroT temp1]];

%% Plot f2f data
FigTrace = figure();
set(FigTrace,'units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
ax1 = subplot('Position',[0.05 0.06 0.90 0.43]); hold on;
for i=1:size(f2fTime,2)
    col = repmat((i-1)*colorStep,1,3);
    plot(f2fTime(:,i),f2fDataS(:,i),'Color',colMap(i,:));
end
grid on;

ax2 = subplot('Position',[0.05 0.54 0.90 0.43]); hold on;
for i=1:size(f2fTime,2)
    col = repmat((i-1)*colorStep,1,3);
    plot(f2fTime(:,i),f2fDataP(:,i),'Color',colMap(i,:));
end
grid on;
linkaxes([ax1,ax2],'x');
xlim([0,0.3e-4]);
title('F2F');

%% Plot Al data
FigTrace = figure();
set(FigTrace,'units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
ax1 = subplot('Position',[0.05 0.06 0.90 0.43]); hold on;
for i=1:size(f2fTime,2)
    col = repmat((i-1)*colorStep,1,3);
    plot(AlTime(:,i),AlDataS(:,i),'Color',colMap(i,:));
end
grid on;

ax2 = subplot('Position',[0.05 0.54 0.90 0.43]); hold on;
for i=1:size(f2fTime,2)
    col = repmat((i-1)*colorStep,1,3);
    plot(AlTime(:,i),AlDataP(:,i),'Color',colMap(i,:));
end
grid on;
linkaxes([ax1,ax2],'x');
xlim([0,0.3e-4]);
title('Al');