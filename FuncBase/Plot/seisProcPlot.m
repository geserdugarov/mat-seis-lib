%% *seisProcPlot*
% Plotting seismograms

%% 
% *Input:* 
% time      - [N]            time array
% Data      - [N, traceNum]  traces
% maxval    - [traceNum]     max values for trace normalizing
% waveNum   - [waveNum]      number of waves
% maxpos    - [traceNum, 3]  first break data
% cross0    - [traceNum, 2]  zero cross data
% pickspath - [string]       path to file with manual pick data
% style     - [ ]            style preferences
%   { 1 - traces style, 2 - automatic pick style, 3 - manual pick style, 
%     4 - zero cross marker style, 
%     5 - x axis label, 6 - y axis label, 7 - title,
%     8 - y axis tick labels,
%     9 - number of traces plotted in window }
% ener      - [traceNum, 1]  estimated wave energy

%%
% *Output:*
% fig    - [id] figure ID

%%
% *Author:* Geser Dugarov 2018

%%

function [fig] = seisProcPlot(time, Data, maxval, waveNum, maxpos, cross0, pickspath, style, ener)

% [thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));

reverseplot = true; % small offset at the top

fig = figure;
set(fig,'units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
ax1 = subplot('Position',[0.05 0.09 0.80 0.88]); hold on;
% plot traces
if reverseplot
    Data = -Data;
end
if maxval > 0
    tempData = repmat(1:size(Data,2),size(Data,1),1) + Data./repmat(maxval,size(Data,1),1)/2;
else
    tempData = repmat(1:size(Data,2),size(Data,1),1) + Data./2;
end
plot(repmat(time,1,size(Data,2)),tempData,style{1});

for waveID = 1:waveNum
    maxposT = maxpos(:,(waveID-1)*3+1:(waveID-1)*3+3);
    cross0T = cross0(:,(waveID-1)*2+1:(waveID-1)*2+2);
    if reverseplot
        maxposT(:,2) = -maxposT(:,2);
    end
    % plot automatic picks
    tempData = (1:size(maxposT,1))' + maxposT(:,2)./maxval'/2;
    plot(maxposT(:,1),tempData,style{2});
    % plot manual picks
    pos = find(maxposT(:,3));
    plot(maxposT(pos,1),tempData(pos,1),style{3});
    % plot zero crosses
    tempData = (1:size(cross0T,1))';
    plot(cross0T(:,1),tempData,style{4}{waveID});
    plot(cross0T(:,2),tempData,style{4}{waveID});
end

% plot settings
xlim([0,max(time)]);
ylim([0,style{9}]);
LinePlotExplorer_polyfit(pickspath);
set(gca,'Ytick',1:size(Data,2)); set(gca,'YtickLabel',style{8});
if reverseplot
    set(gca,'Ydir','reverse');
end
grid on; xlabel(style{5}); ylabel(style{6});
title(style{7});

ax2 = subplot('Position',[0.88 0.09 0.1 0.88]); hold on;
for waveID = 1:waveNum
    enerT = ener(:,waveID);
    tempData = (1:size(ener,1))';
    plot(enerT(:,1),tempData,style{4}{waveID});
end
set(gca,'Ytick',1:size(Data,2)); set(gca,'YtickLabel','');
if reverseplot
    set(gca,'Ydir','reverse');
end
grid on; xlabel('energy');

linkaxes([ax1,ax2],'y');

end