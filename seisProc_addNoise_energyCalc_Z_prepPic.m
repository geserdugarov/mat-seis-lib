%% Tool for processing seismograms
%  add noise
%  find first breaks
%  estimate wave energy
%  Geser Dugarov, 2019

clear variables;
close all;
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
fprintf('    is running... \n');

localMatlabFolder1 = 'FuncBase';  addpath(genpath(localMatlabFolder1));
localMatlabFolder1 = '_VG Private M-Functions';  addpath(genpath(localMatlabFolder1));

% settings for seismogram processing
% 'SP-az000','SP-az030','SP-az060','SP-az090','SP-az120','SP-az150','SP-az180','SP-az210','SP-az240','SP-az270','SP-az300','SP-az330'
%     'ID',      {{'VCH-az000','VCH-az015','VCH-az030','VCH-az045','VCH-az060','VCH-az075','VCH-az090'}}, ...  % seismogram IDs
sett = struct(...
    'ID',      {{'VCH-az000','VCH-az030','VCH-az060'}}, ...  % seismogram IDs
    'waveNum',  1, ...       % number of analysed waves
    'useWind',  true, ... % use window for searching or not
    'searchWind',  [0.4 0.8], ... % time windows for wave searching
    'useFixed',  true, ... % use fixed window for energy estimation or not
    'enerWind',  [0.04,0.04], ... % use fixed window for energy estimation or not
    'marker',      {{'+k','xk','ok','*r'}}, ... % marker settings for different waves
    'firstVal',    0.8, ...     % search first sett.firstVal*max (in window)
    'firstValMod', true, ...    % if true, then max from |x|, else from x
    'maxWind',  0.040, ... % s, window after range*max
    'intWind',  0.008, ... % s, window length for interpolatation and finding max
    'polyfitdegr', 3, ...      % degree of interpolation polynom
    'polyfitdisc', 1.0e-5, ... % interpolation discretization
    'SNR', 2, ...    % signal to noise ratio for noise adding
    'noiseWind',  [0.0 0.3], ... % time window for noise analysis
    'generateID', 1:1, ... % IDs for generated traces
    'noiseFilter', false,                   ... % filter noise or not
    'filtFreqRange', [0 0 60 80] ... % filter before processing ./^^\.
    );

%% Load data for noise generation
noise = load('data\seismic_noise.mat');

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
    foldID = ['VCH-SNR',num2str(sett.SNR),'_',foldID];
    if ~exist([thisFolderName,'\outputs\',foldID],'dir')
        mkdir([thisFolderName,'\outputs\',foldID]);
    end
    fileID = fopen([thisFolderName,'\outputs\',foldID,'\info.txt'],'w');
    % define a number of traces
    seisNum = numel(sett.ID);
    load([thisFolderName,'\data\',sett.ID{1},'-off50-50-2000-PP.mat']);
    traceNum  = size(seisZ,2);
    
    levelSignalInverse = zeros(traceNum*seisNum,sett.waveNum);
    levelSignalReal = zeros(traceNum*seisNum,sett.waveNum);
    waveNormInverse = zeros(3,traceNum*seisNum,sett.waveNum);
    for seisID = 1:seisNum
        load([thisFolderName,'\data\',sett.ID{seisID},'-off50-50-2000-PP.mat']);
        seisZ0 = seisZ;
        fprintf(fileID,'\nAzimuth: ''%s''\n',sett.ID{seisID});

        traceNum  = size(seisZ0,2);
        traceLeng = numel(time);
        % incident angles for normalization
        incAng = rad2deg(acos( waveNorm(3,:) ));
        % normalization for angles less then 7 degrees
        normPos = incAng < 7;

        tracesInWind = 41;   % number of plotted traces in the window

        % estimating of energy before adding noise
        timeBase = zeros(traceNum,sett.waveNum);
        levelSignal = zeros(traceNum,sett.waveNum);
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
            
            % calculating wave energy before adding noise (Z component)
            pickpath = [thisFolderName,'\data\',sett.ID{seisID},'_manualPick_Z.dat'];
            data = seisZ0;
            maxpos = seisProcFindFB(time, data, settTemp, pickpath, []);
            cross0 = seisProcFindCross(time, data, maxpos);
            timeBase(:,waveID) = cross0(:,2) - cross0(:,1);
            levelSignal(:,waveID) = seisProcCalcRMS(time, data, cross0);
            levelSign = mean(levelSignal(:,waveID));
        end

        % cutting noise, calculating level and adding
        noiseCutZ = zeros(traceLeng, traceNum);
        for id = 1:traceNum
            ind = floor(rand*( numel(noise.t)-traceLeng ) + 1);
            noiseCutZ(:,id) = noise.signal(ind:ind+traceLeng-1, 3);
        end
        levelNoise = 3*mean(std(noiseCutZ)); % !!! 99.7% of the values in the range (mu-3*sigma,mu+3*sigma)
        noiseCutZ = noiseCutZ*levelSign/sett.SNR/levelNoise;
        fprintf(fileID,'\n');
        for waveID = 1:sett.waveNum
            fprintf(fileID,'SNR (''%.0f'' wave) = ''%.1f''\n', waveID, sett.SNR);
        end
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
                tempZ = fft(seisZ(:,trace));
                tempZ = real(ifft(tempZ.*filt'));
                seisZ(:,trace) = tempZ*(max(seisZ(:,trace))/max(tempZ));
                clear df freq filt pos1 pos2 pos3 pos4 temp;
            end
        end
        
        % trace normalization with relative amplitude saving (from trace to trace)
        gain = 3;
        maxval = max(max(abs(seisZ)))/gain;
        maxval = repmat(maxval,1,size(seisZ,2));

        % estimating of energy after adding noise
        maxpos = zeros(traceNum,sett.waveNum*3); % for saving windows settings
        cross0 = zeros(traceNum,sett.waveNum*2);
        levelSignalEner = zeros(traceNum,sett.waveNum);
        levelNoiseEner = zeros(traceNum,sett.waveNum);
        levelSignalCorr = zeros(traceNum,sett.waveNum);
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
            if sett.useFixed
                cross0(:,(waveID-1)*2+1) = maxpos(:,(waveID-1)*3+1) - sett.enerWind(waveID);
                cross0(:,(waveID-1)*2+2) = maxpos(:,(waveID-1)*3+1) + sett.enerWind(waveID);
            else
                cross0(:,(waveID-1)*2+1:(waveID-1)*2+2) = seisProcFindCross(time, data, maxpos(:,(waveID-1)*3+1:(waveID-1)*3+3));
            end
            levelSignalEner(:,waveID) = seisProcCalcEner(time, data, cross0(:,(waveID-1)*2+1:(waveID-1)*2+2));

            % estimation and correcting for additional noise energy
            levelNoiseEner(:,waveID) = seisProcCalcEner(time, data, repmat(sett.noiseWind,traceNum,1));
            levelNoiseEner(:,waveID) = levelNoiseEner(:,waveID)/(sett.noiseWind(2)-sett.noiseWind(1)).*(cross0(:,(waveID-1)*2+2) - cross0(:,(waveID-1)*2+1));
            levelSignalEner(:,waveID) = levelSignalEner(:,waveID) - levelNoiseEner(:,waveID);
  
            % energy couldn't be less then 0
            levelSignalEner(levelSignalEner(:,waveID)<0, waveID) = 0;
            
            % convert to relative amplitudes
            levelSignal = sqrt(levelSignalEner);
            levelSignal(:,waveID) = levelSignal(:,waveID)/mean(levelSignal(normPos,waveID));
        end

        switch seisID
            case 1
                seisZ1 = seisZ;
                maxposCorr(:,:,1) = maxpos;
                cross0Corr(:,:,1) = cross0;
            case 2
                seisZ2 = seisZ;
                maxposCorr(:,:,2) = maxpos;
                cross0Corr(:,:,2) = cross0;
            case 3
                seisZ3 = seisZ;
                maxposCorr(:,:,3) = maxpos;
                cross0Corr(:,:,3) = cross0;
        end


    end

end
    % style     - [ ]            style preferences
    %   { 1 - traces style, 2 - automatic pick style, 3 - manual pick style, 
    %     4 - zero cross marker style, 
    %     5 - x axis label, 6 - y axis label, 7 - title,
    %     8 - y axis tick labels,
    %     9 - number of traces plotted in window }
    fig = figure;
    set(fig,'units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
    
    style = {'b','+k','+c',sett.marker,'врем€, с','вынос, м','јзимут 0^{o}', ...
             offset,tracesInWind, [0.3, 1.0]};
    seisProcPlot3CinOne(fig, time, seisZ1, maxval, sett.waveNum, maxposCorr(:,:,1), ...
             cross0Corr(:,:,1), pickpath, style, 1);

    style{1} = 'b'; style{7} = 'јзимут 30^{o}';
    seisProcPlot3CinOne(fig, time, seisZ2, maxval, sett.waveNum, maxposCorr(:,:,2), ...
             cross0Corr(:,:,2), pickpath, style, 2);

    style{1} = 'b'; style{7} = 'јзимут 60^{o}';
    seisProcPlot3CinOne(fig, time, seisZ3, maxval, sett.waveNum, maxposCorr(:,:,3), ...
             cross0Corr(:,:,3), pickpath, style, 3);
    set(fig,'PaperPositionMode','auto');
    saveas(fig,[thisFolderName,'\outputs\seism_example_',sett.ID{seisID}],'png');
        
fprintf('>>> Return from script ''%s'' \n \n', [thisFileName, '.m']);
return