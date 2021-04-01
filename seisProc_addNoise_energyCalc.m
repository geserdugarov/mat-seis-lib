%% Tool for processing seismograms
%  add noise
%  find first breaks
%  estimate wave energy
%  Geser Dugarov, 2018

clear variables;
% close all;
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
    'useFixed',  true, ... % use fixed window for energy estimation or not
    'enerWind',  [0.04,0.04], ... % use fixed window for energy estimation or not
    'marker',      {{'+k','xk','ok'}}, ... % marker settings for different waves
    'firstVal',    0.8, ...     % search first sett.firstVal*max (in window)
    'firstValMod', true, ...    % if true, then max from |x|, else from x
    'maxWind',  0.040, ... % s, window after range*max
    'intWind',  0.008, ... % s, window length for interpolatation and finding max
    'polyfitdegr', 3, ...      % degree of interpolation polynom
    'polyfitdisc', 1.0e-5, ... % interpolation discretization
    'SNR', 500, ...    % signal to noise ratio for noise adding
    'noiseWind',  [0.0 0.4], ... % time window for noise analysis
    'generateID', 1, ... % IDs for generated traces
    'noiseFilter', false,                   ... % filter noise or not
    'filtFreqRange', [0.0 0.0 60 80] ... % filter before processing ./^^\.
    );

%% Load data from azimuth 150 for SNR estimation
ID = 'az150';
load([thisFolderName,'\data\',ID,'_off50-50-2000_PP-PS.mat']);
seisX0 = seisX; seisY0 = seisY; seisZ0 = seisZ;
traceNum  = size(seisX0,2);
traceLeng = numel(time);
enerSignalClear = zeros(traceNum,3,sett.waveNum);
enerSignalClearAllComp = zeros(traceNum,sett.waveNum);
enerSignalClearWave    = zeros(1,sett.waveNum);
enerSignalClearBase    = zeros(traceNum,3,sett.waveNum);
cross0SaveClear = zeros(traceNum,6,sett.waveNum); % !!!
% choosing wave
waveID = 1;
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
% calculating wave energy before adding noise (X component)
pickpath = [thisFolderName,'\data\',ID,'_manualPick_X.dat'];
data = seisX0;
maxpos = seisProcFindFB(time, data, settTemp, pickpath, []);
cross0 = seisProcFindCross(time, data, maxpos);
enerSignalClearBase(:,1,waveID) = cross0(:,2) - cross0(:,1);
enerSignalClear(:,1,waveID) = seisProcCalcEner(time, data, cross0);
% calculating wave energy before adding noise (Y component)
pickpath = [thisFolderName,'\data\',ID,'_manualPick_Y.dat'];
data = seisY0;
maxpos = seisProcFindFB(time, data, settTemp, pickpath, []);
cross0 = seisProcFindCross(time, data, maxpos);
enerSignalClearBase(:,2,waveID) = cross0(:,2) - cross0(:,1);
enerSignalClear(:,2,waveID) = seisProcCalcEner(time, data, cross0);
% calculating wave energy before adding noise (Z component)
pickpath = [thisFolderName,'\data\',ID,'_manualPick_Z.dat'];
data = seisZ0;
maxpos = seisProcFindFB(time, data, settTemp, pickpath, []);
cross0 = seisProcFindCross(time, data, maxpos);
enerSignalClearBase(:,3,waveID) = cross0(:,2) - cross0(:,1);
enerSignalClear(:,3,waveID) = seisProcCalcEner(time, data, cross0);
% estimation of signal level
enerSignalClearAllComp(:,waveID) = enerSignalClear(:,1,waveID) + ...
              enerSignalClear(:,2,waveID) + enerSignalClear(:,3,waveID);
enerSignalClearWave(1,waveID) = mean(enerSignalClearAllComp(:,waveID));
% using first wave as primary wave
levelSign = enerSignalClearWave(1,waveID);

%% Load data for noise generation
noise = load('data\seismic_noise.mat');
% calculating noise level
enerNoise = sum(sum(noise.signal(:,:).^2));
% energy level of noise for the whole trace -> for mean window
levelNoise = mean(enerNoise)/(noise.t(end)-noise.t(1))*mean(mean(enerSignalClearBase(:,:,waveID)));
corrNoise = sqrt(levelSign/sett.SNR/levelNoise);
% correcting noise level
noise.signal = noise.signal.*corrNoise;
enerNoise = sum(sum(noise.signal(:,:).^2));
levelNoise = mean(enerNoise)/(noise.t(end)-noise.t(1))*mean(mean(enerSignalClearBase(:,:,waveID)));
fprintf('SNR (''%.0f'' wave) = ''%.1f''\n', waveID, levelSign/levelNoise);

%% Generating data with noise
for generateID = sett.generateID
    % creating folder name
    if numel(num2str(generateID)) == 1
        foldID = ['00',num2str(generateID)];
    elseif numel(num2str(generateID)) == 2
        foldID = ['0',num2str(generateID)];
    else
        foldID = num2str(generateID);
    end
    % adding noise level in folder name
    foldID = ['noise',num2str(round(100/sett.SNR)),'_',foldID];
    
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
        seisX0 = seisX; seisY0 = seisY; seisZ0 = seisZ;
        fprintf(fileID,'\nAzimuth: ''%s''\n',sett.ID{seisID});

        traceNum  = size(seisX0,2);
        traceLeng = numel(time);

        tracesInWind = 41;   % number of plotted traces in the window

        % estimating of energy before adding noise
        enerSignalClear = zeros(traceNum,3,sett.waveNum);
        enerSignalClearAllComp = zeros(traceNum,sett.waveNum);
        enerSignalClearWave    = zeros(1,sett.waveNum);
        enerSignalClearBase    = zeros(traceNum,3,sett.waveNum);
        cross0SaveClear = zeros(traceNum,6,sett.waveNum); % !!!
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

            % calculating wave energy before adding noise (X component)
            pickpath = [thisFolderName,'\data\',sett.ID{seisID},'_manualPick_X.dat'];
            data = seisX0;
            maxpos = seisProcFindFB(time, data, settTemp, pickpath, []);
            cross0 = seisProcFindCross(time, data, maxpos);
            enerSignalClearBase(:,1,waveID) = cross0(:,2) - cross0(:,1);
            cross0SaveClear(:,1:2,waveID) = cross0; % !!!
            enerSignalClear(:,1,waveID) = seisProcCalcEner(time, data, cross0);
            % calculating wave energy before adding noise (Y component)
            pickpath = [thisFolderName,'\data\',sett.ID{seisID},'_manualPick_Y.dat'];
            data = seisY0;
            maxpos = seisProcFindFB(time, data, settTemp, pickpath, []);
            cross0 = seisProcFindCross(time, data, maxpos);
            enerSignalClearBase(:,2,waveID) = cross0(:,2) - cross0(:,1);
            cross0SaveClear(:,3:4,waveID) = cross0; % !!!
            enerSignalClear(:,2,waveID) = seisProcCalcEner(time, data, cross0);
            % calculating wave energy before adding noise (Z component)
            pickpath = [thisFolderName,'\data\',sett.ID{seisID},'_manualPick_Z.dat'];
            data = seisZ0;
            maxpos = seisProcFindFB(time, data, settTemp, pickpath, []);
            cross0 = seisProcFindCross(time, data, maxpos);
            enerSignalClearBase(:,3,waveID) = cross0(:,2) - cross0(:,1);
            cross0SaveClear(:,5:6,waveID) = cross0; % !!!
            enerSignalClear(:,3,waveID) = seisProcCalcEner(time, data, cross0);

            % estimation of signal level
            enerSignalClearAllComp(:,waveID) = enerSignalClear(:,1,waveID) + ...
                          enerSignalClear(:,2,waveID) + enerSignalClear(:,3,waveID);
            enerSignalClearWave(1,waveID) = mean(enerSignalClearAllComp(:,waveID));
        end

        % cutting noise
        noiseCutX = zeros(traceLeng, traceNum);
        noiseCutY = zeros(traceLeng, traceNum);
        noiseCutZ = zeros(traceLeng, traceNum);
        for id = 1:traceNum
            ind = floor(rand*( numel(noise.t)-traceLeng ) + 1);
            noiseCutX(:,id) = noise.signal(ind:ind+traceLeng-1, 1);
            noiseCutY(:,id) = noise.signal(ind:ind+traceLeng-1, 2);
            noiseCutZ(:,id) = noise.signal(ind:ind+traceLeng-1, 3);
        end

        % calculating noise level
        enerNoise = sum(noiseCutX.^2) + sum(noiseCutY.^2) + sum(noiseCutZ.^2);
        levelNoise = mean(enerNoise)/(time(end)-time(1))*mean(mean(enerSignalClearBase(:,:,1)));
        fprintf(fileID,'\n');
        for waveID = 1:sett.waveNum
            fprintf(fileID,'SNR (''%.0f'' wave) = ''%.1f''\n', waveID, enerSignalClearWave(1,waveID)/levelNoise);
        end

        % adding noise
        seisX = seisX0 + noiseCutX;
        seisY = seisY0 + noiseCutY;
        seisZ = seisZ0 + noiseCutZ;
        
        % filtration
        if sett.noiseFilter
            for trace = 1:traceNum
                df = 1/(time(end)-time(1));
                freq = 0:df:df*(numel(time)-1);
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
                tempX = fft(seisX(:,trace));
                tempY = fft(seisY(:,trace));
                tempZ = fft(seisZ(:,trace));
                tempX = real(ifft(tempX.*filt'));
                tempY = real(ifft(tempY.*filt'));
                tempZ = real(ifft(tempZ.*filt'));
                seisX(:,trace) = tempX*(max(seisX(:,trace))/max(tempX));
                seisY(:,trace) = tempY*(max(seisY(:,trace))/max(tempY));
                seisZ(:,trace) = tempZ*(max(seisZ(:,trace))/max(tempZ));
                clear df freq filt pos1 pos2 pos3 pos4 temp;
            end
        end
        
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
            if sett.useFixed
                cross0(:,(waveID-1)*2+1,1) = maxpos(:,(waveID-1)*3+1,1) - sett.enerWind(waveID);
                cross0(:,(waveID-1)*2+2,1) = maxpos(:,(waveID-1)*3+1,1) + sett.enerWind(waveID);
            else
                cross0(:,(waveID-1)*2+1:(waveID-1)*2+2,1) = seisProcFindCross(time, data, maxpos(:,(waveID-1)*3+1:(waveID-1)*3+3,1));
            end
            enerSignalNoise(:,1,waveID) = seisProcCalcEner(time, data, cross0(:,(waveID-1)*2+1:(waveID-1)*2+2,1));
            % calculating wave energy after adding noise (Y component)
            pickpath = [thisFolderName,'\data\',sett.ID{seisID},'_manualPick_Y.dat'];
            data = seisY;
            maxpos(:,(waveID-1)*3+1:(waveID-1)*3+3,2) = seisProcFindFB(time, data, settTemp, pickpath, []);
            if sett.useFixed
                cross0(:,(waveID-1)*2+1,2) = maxpos(:,(waveID-1)*3+1,2) - sett.enerWind(waveID);
                cross0(:,(waveID-1)*2+2,2) = maxpos(:,(waveID-1)*3+1,2) + sett.enerWind(waveID);
            else
                cross0(:,(waveID-1)*2+1:(waveID-1)*2+2,2) = seisProcFindCross(time, data, maxpos(:,(waveID-1)*3+1:(waveID-1)*3+3,2));
            end
            enerSignalNoise(:,2,waveID) = seisProcCalcEner(time, data, cross0(:,(waveID-1)*2+1:(waveID-1)*2+2,2));
            % calculating wave energy after adding noise (Z component)
            pickpath = [thisFolderName,'\data\',sett.ID{seisID},'_manualPick_Z.dat'];
            data = seisZ;
            maxpos(:,(waveID-1)*3+1:(waveID-1)*3+3,3) = seisProcFindFB(time, data, settTemp, pickpath, []);
            if sett.useFixed
                cross0(:,(waveID-1)*2+1,3) = maxpos(:,(waveID-1)*3+1,3) - sett.enerWind(waveID);
                cross0(:,(waveID-1)*2+2,3) = maxpos(:,(waveID-1)*3+1,3) + sett.enerWind(waveID);
            else
                cross0(:,(waveID-1)*2+1:(waveID-1)*2+2,3) = seisProcFindCross(time, data, maxpos(:,(waveID-1)*3+1:(waveID-1)*3+3,3));
            end
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

            % estimation of noise energy
            enerNoise(:,1,waveID) = seisProcCalcEner(time, noiseCutX, repmat(sett.noiseWind,traceNum,1));
            enerNoise(:,2,waveID) = seisProcCalcEner(time, noiseCutY, repmat(sett.noiseWind,traceNum,1));
            enerNoise(:,3,waveID) = seisProcCalcEner(time, noiseCutZ, repmat(sett.noiseWind,traceNum,1));
            enerNoise(:,1,waveID) = enerNoise(:,1,waveID)/( sett.noiseWind(2)-sett.noiseWind(1) ).*...
                ( cross0Corr(:,(waveID-1)*2+2,1)-cross0Corr(:,(waveID-1)*2+1,1) );
            enerNoise(:,2,waveID) = enerNoise(:,2,waveID)/( sett.noiseWind(2)-sett.noiseWind(1) ).*...
                ( cross0Corr(:,(waveID-1)*2+2,2)-cross0Corr(:,(waveID-1)*2+1,2) );
            enerNoise(:,3,waveID) = enerNoise(:,3,waveID)/( sett.noiseWind(2)-sett.noiseWind(1) ).*...
                ( cross0Corr(:,(waveID-1)*2+2,3)-cross0Corr(:,(waveID-1)*2+1,3) );
            
            % correcting for additional noise energy
            enerSignalNoise(:,1,waveID) = enerSignalNoise(:,1,waveID) - enerNoise(:,1,waveID);
            enerSignalNoise(:,2,waveID) = enerSignalNoise(:,2,waveID) - enerNoise(:,2,waveID);
            enerSignalNoise(:,3,waveID) = enerSignalNoise(:,3,waveID) - enerNoise(:,3,waveID);
            
            % collecting energy data from all 3 components
            enerSignalNoiseAllComp(:,waveID) = enerSignalNoise(:,1,waveID)+enerSignalNoise(:,2,waveID)+enerSignalNoise(:,3,waveID);
            enerSignalNoiseOnlyZ(:,waveID)   = enerSignalNoise(:,3,waveID);
            
            % energy couldn't be less then 0
            enerSignalNoiseAllComp(enerSignalNoiseAllComp(:,waveID)<0, waveID) = 0;

            % error estimation
            errMean(:,waveID) = mean(abs(enerSignalClear(:,:,waveID) - enerSignalNoise(:,:,waveID))) ./ mean(enerSignalClear(:,:,waveID));
            errRMS(:,waveID)  = rms(enerSignalClear(:,:,waveID) - enerSignalNoise(:,:,waveID)) ./ mean(enerSignalClear(:,:,waveID));
            errMeanFinal(1,waveID) = mean(abs(enerSignalClearAllComp(:,waveID) - enerSignalNoiseAllComp(:,waveID))) ./ mean(enerSignalClearAllComp(:,waveID));
            errRMSFinal(1,waveID)  = rms(enerSignalClearAllComp(:,waveID) - enerSignalNoiseAllComp(:,waveID)) ./ mean(enerSignalClearAllComp(:,waveID));

            fprintf(fileID,'\nErrors in energy estimation (comparing noisy data with clear data.)');
            fprintf(fileID,'\nFor wave with ID = ''%d'':', waveID);
            component = ['X','Y','Z'];
            for j = 1:3
                fprintf(fileID,'\n  mean error for ''%s'' component = %.1f%%', component(j), errMean(j,waveID)*100);
            end
            fprintf(fileID,'\n  mean error for all components = %.1f%%\n', errMeanFinal(waveID)*100);
            for j = 1:3
                fprintf(fileID,'\n  RMS error for ''%s'' component = %.1f%%',  component(j), errRMS(j,waveID)*100);
            end
            fprintf(fileID,'\n  RMS error for all components = %.1f%%\n', errRMSFinal(waveID)*100);
        end

        enerSignalClearFinal = enerSignalClearAllComp;
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
            plot(offset,enerSignalClearFinal(:,waveID),'--');
            plot(offset,enerSignalNoiseFinal(:,waveID),style{4}{waveID});
            plot(offset,enerSignalNoiseFinal(:,waveID),style{4}{waveID});
        end
        xlim([0 offset(end)]); xlabel('offset, m');
        ylabel('energy'); title([sett.ID{seisID},', all components']);
        grid on;
        set(fig,'PaperPositionMode','auto');
        saveas(fig,[thisFolderName,'\outputs\',foldID,'\',sett.ID{seisID},'_enerData'],'png');

        save([thisFolderName,'\outputs\',foldID,'\',sett.ID{seisID},'.mat'],'time','seisX','seisY','seisZ','offset','waveNorm','enerSignalClearFinal','enerSignalNoiseFinal');

        for waveID = 1:sett.waveNum
            enerSignalClearInverse((seisID-1)*traceNum+1:seisID*traceNum,waveID)  = enerSignalClearFinal(:,waveID);
            enerSignalNoiseInverse((seisID-1)*traceNum+1:seisID*traceNum,waveID)  = enerSignalNoiseFinal(:,waveID);
            enerSignalNoiseInverseZ((seisID-1)*traceNum+1:seisID*traceNum,waveID) = enerSignalNoiseOnlyZ(:,waveID);
            waveNormInverse(:,(seisID-1)*traceNum+1:seisID*traceNum,waveID) = waveNorm(:,:,waveID);
        end
        
        if numel(sett.generateID) > 2
            close all;
        end
    end
    fclose(fileID);
    
    fig = figure;
    set(fig,'units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
    ax1 = subplot('Position',[0.05 0.09 0.90 0.88]); hold on;
    for waveID = 1:sett.waveNum
        plot(enerSignalClearInverse(:,waveID),'--');
        plot(enerSignalNoiseInverse(:,waveID),style{4}{waveID});
%         plot(enerSignalNoiseInverseZ(:,waveID),'o');
    end
    ylabel('energy'); title(['All azimuths, SNR = ',num2str(sett.SNR)]);
    grid on;
    set(fig,'PaperPositionMode','auto');
    saveas(fig,[thisFolderName,'\outputs\',foldID,'_azim_all_enerData'],'png');
    if numel(sett.generateID) > 2
        close all;
    end
    
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
    fclose(fileID);
end

fprintf('>>> Return from script ''%s'' \n \n', [thisFileName, '.m']);
return