%% Tool for processing seismograms
%  find first breaks
%  estimate wave energy
%  Geser Dugarov, 2018

clear variables;
close all;
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
fprintf('    is running... \n');

localMatlabFolder1 = 'FuncBase';  addpath(genpath(localMatlabFolder1));
localMatlabFolder2 = 'data';  addpath(genpath(localMatlabFolder2));

% settings for seismogram processing
sett = struct(...
    'ID',      {{'az060','az105','az150'}}, ...  % seismogram IDs      ,'az075','az090','az105','az120','az135','az150'        ,'az090','az120','az150' 
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

% define a number of traces
seisNum = numel(sett.ID);
load([sett.ID{1},'_off50-50-2000_PP-PS.mat'],'offset');
traceNum  = numel(offset);
for seisID = 1:seisNum
    load([sett.ID{seisID},'_off50-50-2000_PP-PS.mat']);
    seisX0 = seisX; seisY0 = seisY; seisZ0 = seisZ;
    fprintf('\nAzimuth: ''%s''\n',sett.ID{seisID});

    traceNum  = size(seisX0,2);
    traceLeng = numel(time);

    tracesInWind = traceNum+1;   % number of plotted traces in the window

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
    enerSignalNoiseAllComp = zeros(traceNum, sett.waveNum);
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
        pickpath = [thisFolderName,'\',sett.ID{seisID},'_manualPick_X.dat'];
        data = seisX;
        maxpos(:,(waveID-1)*3+1:(waveID-1)*3+3,1) = seisProcFindFB(time, data, settTemp, pickpath, []);
        cross0(:,(waveID-1)*2+1:(waveID-1)*2+2,1) = seisProcFindCross(time, data, maxpos(:,(waveID-1)*3+1:(waveID-1)*3+3,1));
        enerSignalNoise(:,1,waveID) = seisProcCalcEner(time, data, cross0(:,(waveID-1)*2+1:(waveID-1)*2+2,1));
        % calculating wave energy after adding noise (Y component)
        pickpath = [thisFolderName,'\',sett.ID{seisID},'_manualPick_Y.dat'];
        data = seisY;
        maxpos(:,(waveID-1)*3+1:(waveID-1)*3+3,2) = seisProcFindFB(time, data, settTemp, pickpath, []);
        cross0(:,(waveID-1)*2+1:(waveID-1)*2+2,2) = seisProcFindCross(time, data, maxpos(:,(waveID-1)*3+1:(waveID-1)*3+3,2));
        enerSignalNoise(:,2,waveID) = seisProcCalcEner(time, data, cross0(:,(waveID-1)*2+1:(waveID-1)*2+2,2));
        % calculating wave energy after adding noise (Z component)
        pickpath = [thisFolderName,'\',sett.ID{seisID},'_manualPick_Z.dat'];
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

        % collecting energy data from all 3 components
        enerSignalNoiseAllComp(:,waveID) = enerSignalNoise(:,1,waveID)+enerSignalNoise(:,2,waveID)+enerSignalNoise(:,3,waveID);
    end

    enerSignalNoiseFinal = enerSignalNoiseAllComp;

    % style     - [ ]            style preferences
    %   { 1 - traces style, 2 - automatic pick style, 3 - manual pick style, 
    %     4 - zero cross marker style, 
    %     5 - x axis label, 6 - y axis label, 7 - title,
    %     8 - y axis tick labels,
    %     9 - number of traces plotted in window }
    fig = figure;
    set(fig,'units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
    
    style = {'r','+k','+c',sett.marker,'время, с','вынос, м','X', ...
             offset,tracesInWind, [0.5, 1.1]};
    seisProcPlot3CinOne(fig, time, seisX, maxval, sett.waveNum, maxposCorr(:,:,1), ...
             cross0Corr(:,:,1), pickpath, style, 1);

    style{1} = 'g'; style{7} = 'Y';
    seisProcPlot3CinOne(fig, time, seisY, maxval, sett.waveNum, maxposCorr(:,:,2), ...
             cross0Corr(:,:,2), pickpath, style, 2);

    style{1} = 'b'; style{7} = 'Z';
    seisProcPlot3CinOne(fig, time, seisZ, maxval, sett.waveNum, maxposCorr(:,:,3), ...
             cross0Corr(:,:,3), pickpath, style, 3);
    set(fig,'PaperPositionMode','auto');
    saveas(fig,[thisFolderName,'\outputs\seism_example_',sett.ID{seisID}],'png');

end

fprintf('>>> Return from script ''%s'' \n \n', [thisFileName, '.m']);
return