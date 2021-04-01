%% Tool for plotting text files with ultrasonic experiment data
%  Geser Dugarov, 2018
%% Set paths and directories
clear variables; close all;
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
fprintf('    is running... \n');
addpath(genpath('FuncBase'));

%% Settings
% path to the folder with data
% path = 'c:\Geser\hydrates_data\20181120 - тестовые измерения в ИЯФ с заморозкой\';
% ID = '20181120'; plotTraceNum = 20;
% path = 'd:\Dropbox\Projects\проект 2018-2019 СО РАН интеграционный (газогидраты)\томография в ИЯФ\компактная акустика\2019.10.13-15 - съемка в ИЯФ\20191019\';
% ID = '20191019'; plotTraceNum = 100;
% path = 'd:\Dropbox\Projects\проект 2018-2019 СО РАН интеграционный (газогидраты)\томография в ИЯФ\компактная акустика\2019.10.13-15 - съемка в ИЯФ\20191015\';
% ID = '20191015'; plotTraceNum = 40;
% path = 'd:\Dropbox\Projects\проект 2018-2019 СО РАН интеграционный (газогидраты)\томография в ИЯФ\компактная акустика\2019.10.13-15 - съемка в ИЯФ\20191018\';
% ID = '20191018'; plotTraceNum = 100;
path = '\\ipgg\project\Проект 14-17\2020\202002 Февраль\20200205\';
ID = '20200207'; plotTraceNum = 100;

gain = 10; % amplification for traces plotting

numLines = 10053; % number of lines in text files
skipLines = 3;   % number of lines that should be skipped
step = 1;        % step for skipping files (1 - without skip)

F = dir([path,'*_PS_*.txt']); %Выбор данных по идентификатору
fileID = 1:step:length(F);
tracenum = numel(fileID);
valuenum = numLines-skipLines;
Time  = nan(valuenum,tracenum);
DataP = nan(valuenum,tracenum);
DataS = nan(valuenum,tracenum);
count = 1;
for i = fileID
    fname = F(i).name;
    temp = importdata([path,fname],'\t',skipLines);
    temp = temp.data;
    Time(:,count)  = temp(:,1);  
    DataP(:,count) = temp(:,3);
    DataS(:,count) = temp(:,6);  
    count = count+1;
end

%% Plotting traces
Data = DataP;
FigTrace = figure();
set(FigTrace,'units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
subplot('Position',[0.05 0.06 0.90 0.88]); hold on;
% max values for normalization of traces
maxval = repmat(max(max(abs(Data))),1,size(Data,2));
tempData = repmat(1:size(Data,2),size(Data,1),1) + gain*Data./repmat(maxval,size(Data,1),1)/2;
plot(Time,tempData,'k');
clear maxval tempData pos;
xlim([0,Time(end,1)]);
ylim([0,plotTraceNum]);
pickspath = [path,'procTextData',ID,'.txt'];
LinePlotExplorer_polyfit(pickspath);
% figure settings for traces
set(gca,'Ytick',1:tracenum); grid on; xlabel('times, s'); ylabel('trace number'); title([ID,', P wave'],'Interpreter','none');
clear pickspath;

Data = DataS;
FigTrace = figure();
set(FigTrace,'units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
subplot('Position',[0.05 0.06 0.90 0.88]); hold on;
% max values for normalization of traces
maxval = repmat(max(max(abs(Data))),1,size(Data,2));
tempData = repmat(1:size(Data,2),size(Data,1),1) + gain*Data./repmat(maxval,size(Data,1),1)/2;
plot(Time,tempData,'k');
clear maxval tempData pos;
xlim([0,Time(end,1)]);
ylim([0,plotTraceNum]);
pickspath = [path,'procTextData',ID,'.txt'];
LinePlotExplorer_polyfit(pickspath);
% figure settings for traces
set(gca,'Ytick',1:tracenum); grid on; xlabel('times, s'); ylabel('trace number'); title([ID,', S wave'],'Interpreter','none');
clear pickspath;


fprintf('\n>>> Return from script ''%s'' \n \n', [thisFileName, '.m']);
return;