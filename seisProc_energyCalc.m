%% Tool for processing seismograms
%  add noise
%  find first breaks
%  estimate wave energy
%  Geser Dugarov, 2018

clear variables;
close all;
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
fprintf('    is running... \n');

localMatlabFolder1 = 'FuncBase';  addpath(genpath(localMatlabFolder1));
localMatlabFolder1 = '_VG Private M-Functions';  addpath(genpath(localMatlabFolder1));

% settings for seismogram processing
sett = struct(...
    'ID',      {{'az060'}}, ...  % seismogram IDs      ,'az075','az090','az105','az120','az135','az150'        ,'az090','az120','az150'
    'waveNum',  2, ...       % number of analysed waves
    'useWind',  true, ... % use window for searching or not
    'searchWind',  [0.5 0.8; 0.8 1.2], ... % time windows for wave searching
    'marker',      {{'+k','xk','ok'}}, ... % marker settings for different waves
    'firstVal',    0.8, ...     % search first sett.firstVal*max (in window)
    'firstValMod', true, ...    % if true, then max from |x|, else from x
    'maxWind',  0.040, ... % s, window after range*max
    'intWind',  0.008, ... % s, window length for interpolatation and finding max
    'polyfitdegr', 3, ...      % degree of interpolation polynom
    'polyfitdisc', 1.0e-5 ... % interpolation discretization
    );

%% Processing data
foldID = 'energyCalcOnly';
if ~exist([thisFolderName,'\outputs\',foldID],'dir')
    mkdir([thisFolderName,'\outputs\',foldID]);
end
fileID = fopen([thisFolderName,'\outputs\',foldID,'\info.txt'],'w');
% define a number of traces
seisNum = numel(sett.ID);
load([thisFolderName,'\data\',sett.ID{1},'_off50-50-2000_PP-PS.mat'],'seisX');
traceNum  = size(seisX,2);
enerSignalClearInverse  = zeros(traceNum*seisNum,sett.waveNum);
enerSignalNoiseInverse  = zeros(traceNum*seisNum,sett.waveNum);
enerSignalNoiseInverseZ = zeros(traceNum*seisNum,sett.waveNum);
waveNormInverse = zeros(3,traceNum*seisNum,sett.waveNum);
for seisID = 1:seisNum
    load([thisFolderName,'\data\',sett.ID{seisID},'_off50-50-2000_PP-PS.mat']);
    fprintf(fileID,'\nAzimuth: ''%s''\n',sett.ID{seisID});

    traceNum  = size(seisX,2);
    traceLeng = numel(time);
    tracesInWind = 41;   % number of plotted traces in the window
    % trace normalization with relative amplitude saving
    % (from trace to trace, from component to component)
    gain = 5;
    maxval = max(max(abs([seisX seisY seisZ])))/gain;
    maxval = repmat(maxval,1,size(seisX,2));

    % estimating of energy after adding noise
    maxpos    = zeros(traceNum,sett.waveNum*3,3); % for saving windows settings
    cross0    = zeros(traceNum,sett.waveNum*2,3);
    maxEnerPos = zeros(1,sett.waveNum);
    maxposCorr = zeros(traceNum,sett.waveNum*3,3);
    enerSignalNoise = zeros(traceNum,3,sett.waveNum); % for saving energy data with noise
    enerNoise = zeros(traceNum,3,sett.waveNum); % for saving energy of noise
    enerSignalNoiseAllComp = zeros(traceNum, sett.waveNum);
    errMean      = zeros(3,sett.waveNum); % for error estimation
    errRMS       = zeros(3,sett.waveNum);
    errMeanFinal = zeros(1,sett.waveNum);
    errRMSFinal  = zeros(1,sett.waveNum);
    for waveID = 1:sett.waveNum
        % creating temporary structure with settings
        settTemp = struct(...
            'useWind',  true, ...
            'searchWind',  sett.searchWind(waveID,:), ...
            'firstVal',    sett.firstVal, ...
            'firstValMod', sett.firstValMod, ...
            'maxWind',  sett.maxWind, ...
            'intWind',  sett.intWind, ...
            'polyfitdegr', sett.polyfitdegr, ...
            'polyfitdisc', sett.polyfitdisc ...
            );

        % calculating wave energy after adding noise (X component)
        pickpath = [thisFolderName,'\data\',sett.ID{seisID},'_manualPick_X.dat'];
        data = seisX;
        maxpos(:,(waveID-1)*3+1:(waveID-1)*3+3,1) = seisProcFindFB(time, data, settTemp, pickpath, []);
        cross0(:,(waveID-1)*2+1:(waveID-1)*2+2,1) = seisProcFindCross(time, data, maxpos(:,(waveID-1)*3+1:(waveID-1)*3+3,1));
        enerSignalNoise(:,1,waveID) = seisProcCalcEner(time, data, cross0(:,(waveID-1)*2+1:(waveID-1)*2+2,1));
        % calculating wave energy after adding noise (Y component)
        pickpath = [thisFolderName,'\data\',sett.ID{seisID},'_manualPick_Y.dat'];
        data = seisY;
        maxpos(:,(waveID-1)*3+1:(waveID-1)*3+3,2) = seisProcFindFB(time, data, settTemp, pickpath, []);
        cross0(:,(waveID-1)*2+1:(waveID-1)*2+2,2) = seisProcFindCross(time, data, maxpos(:,(waveID-1)*3+1:(waveID-1)*3+3,2));
        enerSignalNoise(:,2,waveID) = seisProcCalcEner(time, data, cross0(:,(waveID-1)*2+1:(waveID-1)*2+2,2));
        % calculating wave energy after adding noise (Z component)
        pickpath = [thisFolderName,'\data\',sett.ID{seisID},'_manualPick_Z.dat'];
        data = seisZ;
        maxpos(:,(waveID-1)*3+1:(waveID-1)*3+3,3) = seisProcFindFB(time, data, settTemp, pickpath, []);
        cross0(:,(waveID-1)*2+1:(waveID-1)*2+2,3) = seisProcFindCross(time, data, maxpos(:,(waveID-1)*3+1:(waveID-1)*3+3,3));
        enerSignalNoise(:,3,waveID) = seisProcCalcEner(time, data, cross0(:,(waveID-1)*2+1:(waveID-1)*2+2,3));

        % defining component with max rms amplitude
        [~,maxEnerPos(1,waveID)] = max(mean(enerSignalNoise(:,:,waveID)));

        % saving result with max energy
        temp = maxpos(:,(waveID-1)*3+1,maxEnerPos(1,waveID));
        temp = [(1:numel(temp))' temp];
        cross0save = cross0(:,(waveID-1)*2+1:(waveID-1)*2+2,maxEnerPos(1,waveID));

        % changing zero cross data for all components according max energy
        data = seisX;
        maxposCorr(:,(waveID-1)*3+1:(waveID-1)*3+3,1) = seisProcFindFB(time, data, settTemp, pickpath, temp);
        cross0Corr(:,(waveID-1)*2+1:(waveID-1)*2+2,1) = cross0(:,(waveID-1)*2+1:(waveID-1)*2+2,maxEnerPos(1,waveID));
        enerSignalNoise(:,1,waveID) = seisProcCalcEner(time, data, cross0Corr(:,(waveID-1)*2+1:(waveID-1)*2+2,1));
        data = seisY;
        maxposCorr(:,(waveID-1)*3+1:(waveID-1)*3+3,2) = seisProcFindFB(time, data, settTemp, pickpath, temp);
        cross0Corr(:,(waveID-1)*2+1:(waveID-1)*2+2,2) = cross0(:,(waveID-1)*2+1:(waveID-1)*2+2,maxEnerPos(1,waveID));
        enerSignalNoise(:,2,waveID) = seisProcCalcEner(time, data, cross0Corr(:,(waveID-1)*2+1:(waveID-1)*2+2,2));
        data = seisZ;
        maxposCorr(:,(waveID-1)*3+1:(waveID-1)*3+3,3) = seisProcFindFB(time, data, settTemp, pickpath, temp);
        cross0Corr(:,(waveID-1)*2+1:(waveID-1)*2+2,3) = cross0(:,(waveID-1)*2+1:(waveID-1)*2+2,maxEnerPos(1,waveID));
        enerSignalNoise(:,3,waveID) = seisProcCalcEner(time, data, cross0Corr(:,(waveID-1)*2+1:(waveID-1)*2+2,3));

        cross0SaveNoise = cross0Corr; % !!!
        % interpolation of zero cross data for smooth changing from trace to trace
%         for j = 1:2
%             intpoly = polyfit((1:size(cross0save,1))',cross0save(:,j),sett.polyfitdegr);
%             saveArr = polyval(intpoly,(1:size(cross0save,1))');
%             for k = 1:3
%                 cross0Corr(:,(waveID-1)*2+j,k) = saveArr;
%             end
%         end

        % collecting energy data from all 3 components
        enerSignalNoiseAllComp(:,waveID) = enerSignalNoise(:,1,waveID)+enerSignalNoise(:,2,waveID)+enerSignalNoise(:,3,waveID);
        enerSignalNoiseOnlyZ(:,waveID)   = enerSignalNoise(:,3,waveID);

        % energy couldn't be less then 0
        enerSignalNoiseAllComp(enerSignalNoiseAllComp(:,waveID)<0, waveID) = 0;
    end

    enerSignalNoiseFinal = enerSignalNoiseAllComp;
%     enerSignalNoiseFinal(:,1) = enerSignalNoise(:,3,1);
%     enerSignalNoiseFinal(:,2) = enerSignalNoiseAllComp(:,2);

    % style     - [ ]            style preferences
    %   { 1 - traces style, 2 - automatic pick style, 3 - manual pick style, 
    %     4 - zero cross marker style, 
    %     5 - x axis label, 6 - y axis label, 7 - title,
    %     8 - y axis tick labels,
    %     9 - number of traces plotted in window }
    style = {'r','+k','+c',sett.marker,'times, s','offset, m',[sett.ID{seisID},', X component'], ...
             offset,tracesInWind};
    fig = seisProcPlot(time, seisX, maxval, sett.waveNum, maxposCorr(:,:,1), ...
             cross0Corr(:,:,1), pickpath, style, enerSignalNoise(:,1,:));
    set(fig,'PaperPositionMode','auto');
    saveas(fig,[thisFolderName,'\outputs\',foldID,'\',sett.ID{seisID},'_X'],'png');

    style{1} = 'g'; style{7} = [sett.ID{seisID},', Y component'];
    fig = seisProcPlot(time, seisY, maxval, sett.waveNum, maxposCorr(:,:,2), ...
             cross0Corr(:,:,2), pickpath, style, enerSignalNoise(:,2,:));
    set(fig,'PaperPositionMode','auto');
    saveas(fig,[thisFolderName,'\outputs\',foldID,'\',sett.ID{seisID},'_Y'],'png');

    style{1} = 'b'; style{7} = [sett.ID{seisID},', Z component'];
    fig = seisProcPlot(time, seisZ, maxval, sett.waveNum, maxposCorr(:,:,3), ...
             cross0Corr(:,:,3), pickpath, style, enerSignalNoise(:,3,:));
    set(fig,'PaperPositionMode','auto');
    saveas(fig,[thisFolderName,'\outputs\',foldID,'\',sett.ID{seisID},'_Z'],'png');

    fig = figure;
    set(fig,'units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
    subplot('Position',[0.05 0.09 0.90 0.88]); hold on;
    for waveID = 1:sett.waveNum
        plot(offset,enerSignalNoiseFinal(:,waveID),style{4}{waveID});
        plot(offset,enerSignalNoiseFinal(:,waveID),style{4}{waveID});
    end
    xlim([0 offset(end)]); xlabel('offset, m');
    ylabel('energy'); title([sett.ID{seisID},', all components']);
    grid on;
    set(fig,'PaperPositionMode','auto');
    saveas(fig,[thisFolderName,'\outputs\',foldID,'\',sett.ID{seisID},'_enerData'],'png');

%     save([thisFolderName,'\outputs\',foldID,'\',sett.ID{seisID},'.mat'],'time','seisX','seisY','seisZ','offset','waveNorm','enerSignalNoiseFinal');

    for waveID = 1:sett.waveNum
        enerSignalNoiseInverse((seisID-1)*traceNum+1:seisID*traceNum,waveID)  = enerSignalNoiseFinal(:,waveID);
        enerSignalNoiseInverseZ((seisID-1)*traceNum+1:seisID*traceNum,waveID) = enerSignalNoiseOnlyZ(:,waveID);
        waveNormInverse(:,(seisID-1)*traceNum+1:seisID*traceNum,waveID) = waveNorm(:,:,waveID);
    end

    if numel(sett.ID) > 2
        close all;
    end
end
fclose(fileID);

%% prepare text file for inverse problem
% projection on XY plane
temp = [waveNormInverse(1:2,:,waveID); zeros(1,size(waveNormInverse,2))];
temp = temp./vecnorm(temp);
if ~exist('Outputs', 'dir')
   mkdir('Outputs');
end
fileID = fopen(['Outputs\',foldID,'_forInverseData_PP_PS.txt'],'w');
fprintf(fileID,'%.0f %.4f %.6f %.6f %.4f %.6f %.6f\n', ...
    [rad2deg(acos( temp(1,:) )); rad2deg(acos( waveNormInverse(3,:,1) )); enerSignalClearInverse(:,1)'; enerSignalNoiseInverse(:,1)'; ...
                                 rad2deg(acos( waveNormInverse(3,:,2) )); enerSignalClearInverse(:,2)'; enerSignalNoiseInverse(:,2)']);
fclose('all');

fprintf('>>> Return from script ''%s'' \n \n', [thisFileName, '.m']);
return