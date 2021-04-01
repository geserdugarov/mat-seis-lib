%% *seisProcPlot3CinOne*
% Plotting seismograms

%% 
% *Input:* 
% fig       - object         figure ID
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
%     9 - number of traces plotted in window,
%    10 - xlim parameters}
% posPlot   - [1 2 or 3]     position in figure

%%
% *Output:*
% fig    - [id] figure ID

%%
% *Author:* Geser Dugarov 2018

%%

function [] = seisProcPlot3CinOne(fig, time, Data, maxval, waveNum, maxpos, cross0, pickspath, style, posPlot)

reverseplot = true; % small offset at the top

figure(fig);
% ax1 = subplot('Position',[0.32*posPlot-0.27 0.08 0.27 0.88]); hold on;
ax1 = subplot('Position',[0.32*posPlot-0.27+0.02 0.09 0.27 0.86]); hold on;
Data(:,2:2:end) = []; % !!!
maxval(2:2:end) = []; % !!!
maxpos(2:2:end,:) = []; % !!!
cross0(2:2:end,:) = []; % !!!
% plot traces
if reverseplot
    Data = -Data;
end
if maxval > 0
    tempData = repmat(1:size(Data,2),size(Data,1),1) + Data./repmat(maxval,size(Data,1),1)/2;
else
    tempData = repmat(1:size(Data,2),size(Data,1),1) + Data./2;
end
plot(repmat(time,1,size(Data,2)),tempData,style{1},'LineWidth',2);

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
    plot(maxposT(pos,1),tempData(pos,1),style{3},'LineWidth',2,'MarkerSize',10);
    % plot zero crosses
    tempData = (1:size(cross0T,1))';
    plot(cross0T(:,1),tempData,style{4}{waveID},'LineWidth',2,'MarkerSize',10);
    plot(cross0T(:,2),tempData,style{4}{waveID},'LineWidth',2,'MarkerSize',10);
end

% plot settings
xlim(style{10});
ylim([0,style{9}]);
LinePlotExplorer_polyfit(pickspath);
% set(gca,'Ytick',1:size(Data,2));
% set(gca,'YtickLabel',style{8});
set(gca,'Ytick',1:2:size(Data,2)); % !!!
temp = style{8}; temp(2:2:end) = ''; temp(2:2:end) = ''; % !!!
set(gca,'YtickLabel',temp);
if reverseplot
    set(gca,'Ydir','reverse');
end
grid on;
xlabel(style{5});
if posPlot==1
	ylabel(style{6});
end
set(gca,'FontSize',16); % !!!
ylim([0,21]); % !!!
title(style{7});

end