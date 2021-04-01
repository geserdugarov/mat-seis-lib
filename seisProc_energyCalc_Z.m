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
    'ID',      {{'SP-az000','SP-az030'}}, ...  % seismogram IDs
    'waveNum',  1, ...       % number of analysed waves
    'useWind',  true, ... % use window for searching or not
    'searchWind',  [2.0 3.0], ... % time windows for wave searching
    'marker',      {{'+k','xk','ok','*r'}}, ... % marker settings for different waves
    'firstVal',    0.8, ...     % search first sett.firstVal*max (in window)
    'firstValMod', true, ...    % if true, then max from |x|, else from x
    'maxWind',  0.040, ... % s, window after range*max
    'intWind',  0.008, ... % s, window length for interpolatation and finding max
    'polyfitdegr', 3, ...      % degree of interpolation polynom
    'polyfitdisc', 1.0e-5 ... % interpolation discretization
    );

%% Processing data
foldID = 'SP';
if ~exist([thisFolderName,'\outputs\',foldID],'dir')
    mkdir([thisFolderName,'\outputs\',foldID]);
end
fileID = fopen([thisFolderName,'\outputs\',foldID,'\info.txt'],'w');
% define a number of traces
seisNum = numel(sett.ID);
load([thisFolderName,'\data\',sett.ID{1},'-off100-100-2000-PP.mat'],'seisZ');
traceNum  = size(seisZ,2);
levelSignalInverse  = zeros(traceNum*seisNum,sett.waveNum);
waveNormInverse = zeros(3,traceNum*seisNum,sett.waveNum);
for seisID = 1:seisNum
    load([thisFolderName,'\data\',sett.ID{seisID},'-off100-100-2000-PP.mat']);
    fprintf(fileID,'\nAzimuth: ''%s''\n',sett.ID{seisID});

    traceNum  = size(seisZ,2);
    traceLeng = numel(time);
    tracesInWind = 21;   % number of plotted traces in the window
    % trace normalization with relative amplitude saving
    % (from trace to trace, from component to component)
    gain = 3;
    maxval = max(max(abs(seisZ)))/gain;
    maxval = repmat(maxval,1,size(seisZ,2));

    % estimating of energy after adding noise
    maxpos    = zeros(traceNum,sett.waveNum*3); % for saving windows settings
    cross0    = zeros(traceNum,sett.waveNum*2);
    levelSignalRMS = zeros(traceNum,sett.waveNum); % for saving energy data with noise
    levelSignalRMSCorr = zeros(traceNum,sett.waveNum);
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

        % calculating wave energy after adding noise (Z component)
        pickpath = [thisFolderName,'\data\',sett.ID{seisID},'_manualPick_Z.dat'];
        data = seisZ;
        maxpos(:,(waveID-1)*3+1:(waveID-1)*3+3) = seisProcFindFB(time, data, settTemp, pickpath, []);
        cross0(:,(waveID-1)*2+1:(waveID-1)*2+2) = seisProcFindCross(time, data, maxpos(:,(waveID-1)*3+1:(waveID-1)*3+3));
        levelSignalRMS(:,waveID) = seisProcCalcRMS(time, data, cross0(:,(waveID-1)*2+1:(waveID-1)*2+2));
        % energy couldn't be less then 0
        levelSignalRMS(levelSignalRMS(:,waveID)<0, waveID) = 0;
    end

    % style     - [ ]            style preferences
    %   { 1 - traces style, 2 - automatic pick style, 3 - manual pick style, 
    %     4 - zero cross marker style, 
    %     5 - x axis label, 6 - y axis label, 7 - title,
    %     8 - y axis tick labels,
    %     9 - number of traces plotted in window }
    style = {'b','+k','+c',sett.marker,'times, s','offset, m',[sett.ID{seisID},', Z component'], ...
             offset,tracesInWind};
    fig = seisProcPlot(time, seisZ, maxval, sett.waveNum, maxpos(:,:), ...
             cross0(:,:), pickpath, style, levelSignalRMS(:,:));
    set(fig,'PaperPositionMode','auto');
    saveas(fig,[thisFolderName,'\outputs\',foldID,'\',sett.ID{seisID},'_seisZ'],'png');

    fig = figure;
    set(fig,'units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
    subplot('Position',[0.05 0.09 0.90 0.88]); hold on;
    for waveID = 1:sett.waveNum
        output = levelSignalRMS(:,waveID); output = output/output(1);
        plot(offset,output,style{4}{1});
        output = output./(-waveOut(3,:))'; output = output/output(1);
        plot(offset,output,style{4}{2});
        output = output./(upLayerRefr)'; output = output/output(1);
        levelSignalRMSCorr(:,waveID) = output;
        plot(offset,output,style{4}{3});
        plot(offset,realRefl/realRefl(1),style{4}{3});
    end
    xlim([0 offset(end)]); xlabel('offset, m');
    ylabel('signal level'); title([sett.ID{seisID},', all components']);
    grid on;
    legend('signal level from traces','corrected direction of receiving',...
        'corrected upper layers refractions','real values');
    set(fig,'PaperPositionMode','auto');
    saveas(fig,[thisFolderName,'\outputs\',foldID,'\',sett.ID{seisID},'_dataProc'],'png');

    for waveID = 1:sett.waveNum
        levelSignalInverse((seisID-1)*traceNum+1:seisID*traceNum,waveID) = levelSignalRMSCorr(:,waveID);
        waveNormInverse(:,(seisID-1)*traceNum+1:seisID*traceNum,waveID)  = waveNorm(:,:,waveID);
    end

    if numel(sett.ID) > 2
        close all;
    end
end
fclose(fileID);

%% prepare text file for inverse problem
% projection on XY plane
temp = [waveNormInverse(1:2,:); zeros(1,size(waveNormInverse,2))];
temp = temp./vecnorm(temp);
if ~exist('Outputs', 'dir')
   mkdir('Outputs');
end
fileID = fopen(['Outputs\',foldID,'_forInverseData_PP.txt'],'w');
fprintf(fileID,'%.0f %.2f %.5f\n', ...
    [rad2deg(acos( temp(1,:) )); rad2deg(acos( waveNormInverse(3,:) )); levelSignalInverse(:,1)']);
fclose('all');

fprintf('>>> Return from script ''%s'' \n \n', [thisFileName, '.m']);
return