%% Tool for processing ultrasonic experiment data
%  Velocity and attenuation estimation
%  Max-phase signal is expected
%  Geser Dugarov, 2016

%% Set paths and directories
clear variables; close all;
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
fprintf('    is running... \n');
addpath(genpath('FuncBase'));

%% Data director path
% datadir = '\\ipgg\project\Проект 14-17\Библиотека экспериментов 2016\Данные\'; % path to directory with data
datadir = 'c:\Geser\hydrates_data\'; % path to directory with data

%% Settings for V estimation
settV = struct(...
    'firstVal', 0.35, ...    % search first settV.firstVal*max
    'maxWind',  2.0e-6, ... % s, window after range*max
    'intWind',  0.5e-6, ... % s, window length for interpolatation and finding max
    'polyfitdegr', 3, ...      % degree of interpolation polynom
    'polyfitdisc', 1.0e-9, ... % interpolation discretization
    'noiseFilter', true,                    ... % filter noise or not
    'filtFreqRange', [0.0 0.0 1.5e6 1.6e6]  ... % filter before processing ./^^\.
    );

%% Settings for Q estimation
settQ = struct( ...
    'Qestim', false, ... % true - calculate 1/Q
    'dtw', 0.5e-6, ... % s, the length of window for signal smoothing %0.5e-6
    'freqRangeP', [0.15e6, 0.55e6], ... % Hz, frequency range for P wave %[0.2e6, 0.6e6]
    'freqRangeS', [0.15e6, 0.55e6], ... % Hz, frequency range for S wave %[0.2e6, 0.6e6]
... % timewind, s, [shift to start from max, max window from shifted value] 
... % timewind is not required parameter
... 'timewindP', [1.0e-6, 4.0e-6], ...
... 'timewindS', [5.0e-6, 15.0e-6], ...
    'numOf0CrossP', 2, ... % number of zero crossing after max position for P %3
    'numOf0CrossS', 3, ... % number of zero crossing after max position for S %4
    'plotFreqMax', 1.4e6, ... % Hz, maximum frequency for plotting
    'plotSkipNum', 3 ... % number of traces thar should be skipped
                     ... % while plotting signals, spectrum, Ln(Sp)
	);

%% Settings for saving and plotting
settSavePlot = struct( ...
    'plotTracesNorm', false, ... % true - plot normalized traces, false - real amplitudes
    'plotTraces',   1, ... % 0 or 1 only! 1 - plot traces (used in figID)
    'plotSigAndSp', 0, ... % 0 or 1 only! 1 - plot cutted signals and amplitude spectrums
    'plotLnSp',     0, ... % 0 or 1 only! 1 - plot logarithm of spectrum ratio
    'plotReprOnly', true, ... % plot signals, spectrum and LnSp only for representative traces
    'saveFigsMAT', false, ... % true - save figure mat-files, false - without saving
    'saveFigsPNG', false, ... % true - save figure png-files, false - without saving
    'saveXLS',  false,  ... % true - save xls file, false - without saving
    'plotXLS',  true,  ... % true - plot in xls file, false - without plotting
    'saveSmpForXLS', true ... % true - save only result png-file near xls-file
    );

%% List of sample codes with numbers of bad traces
% 'identificator' [P-wave bad traces] [S-wave bad traces]
% first trace for etalon, second - for aluminum sample
identifList = {
%     'Sthfa' [521] [3:100 521];
%     'Sthfb' [] [3:40 470:481 27:30];
%     'Sthfc' [] [3:70 460:472 19:56 321:323];
%     'Sthfd' [] [3:42 27];
%     'Sthfe' [] [3:43 209:496];
%     'Sthff' [] [3:35 200:482];    
%     'Sthfg' [] [3:352 509:1755 1782:1877 1967:2228];
%     'Sthfh' [] [3:107 740:1526];
%     'Sthfi' [682:683] [3:108 563:683];
%     'Sthfj' [] [];
%     'Sthfk' [] [];
%     'Sthfl' [] [1030:1080];
%     'Sthfm' [] [];
%     'Sthfn' [] [];
%     'Sthfo' [] [3:67 432:469];
%     'Sc01'  [32 67 94 155 190 217 252] [26:33 61:68 88:95 123:129 149:156 609:884 1303:1353 1681:1851];
    'cups' [] [];
    };

traceList = {
%     'Sthfa' [58 343 469];
%     'Sthfb' [5 334 457];
%     'Sthfc' [8 131 166];
%     'Sthfd' [7 93 155];
%     'Sthfe' [490 110 0];
%     'Sthff' [10 0 166];
%     'Sthfg' [1000 0 460];
%     'Sthfh' [10 524 719];
%     'Sthfi' [9 181 0];
%     'Sthfj' [48 170 0];
%     'Sthfk' [650 375 920];
%     'Sthfl' [20 600 0];
%     'Sthfm' [30 410 619];
%     'Sthfn' [50 229 0];
%     'Sthfo' [10 200 0]; %[10 149 200]
%     'Sc01'  [5 0 40];
    'cups' [0 0 0];
    };

%% Saving parameter options in file
if ~exist('outputs','dir')
    mkdir('outputs');
end
save([thisFolderName,'\outputs\settV.mat'],'settV');
save([thisFolderName,'\outputs\settQ.mat'],'settQ');
save([thisFolderName,'\outputs\identifList.mat'],'identifList');
if (settSavePlot.saveFigsMAT || settSavePlot.saveFigsPNG) && ~exist('outputs','dir')
	mkdir('outputs');
end

%% Load data for aluminum
load([thisFolderName,'\data\dataAl.mat'],'dataAl');

if settQ.Qestim
    figID = struct( ...
        'info',     1, ...
        'tracesP',  1+settSavePlot.plotTraces, ...
        'tracesS',  1+2*settSavePlot.plotTraces, ...
        'SigAndSp', 1+2*settSavePlot.plotTraces+settSavePlot.plotSigAndSp, ...
        'lnSp',     1+2*settSavePlot.plotTraces+settSavePlot.plotSigAndSp+settSavePlot.plotLnSp, ...
        'results',  2+2*settSavePlot.plotTraces+settSavePlot.plotSigAndSp+settSavePlot.plotLnSp, ...
        'resVQfT',  3+2*settSavePlot.plotTraces+settSavePlot.plotSigAndSp+settSavePlot.plotLnSp  ...
        );
else
    figID = struct( ...
        'info',     1, ...
        'tracesP',  1+settSavePlot.plotTraces, ...
        'tracesS',  1+2*settSavePlot.plotTraces, ...
        'results',  2+2*settSavePlot.plotTraces, ...
        'resVQfT',  3+2*settSavePlot.plotTraces  ...
        );
end
if size(identifList,1) > 3
    fprintf('WARNING: The number of samples > 3. Figures would be automatically closed after each sample processing.\n');
end
filenamelist = dir(datadir);
for identif = 1:size(identifList,1)
    %% Search and load data
    sampleID = identifList(identif,1); sampleID = sampleID{1};
    fprintf('Processing sample: %s\n', sampleID);
    clear sampname etalname;
    for i = 1:size(filenamelist,1)
        name = filenamelist(i).name;
        searchres1 = strfind(name,[sampleID,'_']);
        searchres2 = strfind(name,[sampleID,'.']);
        searchres3 = strfind(name,[sampleID,' ']);
        searchres = [searchres1 searchres2 searchres3];
        clear searchres1 searchres2 searchres3;
        if ~isempty(searchres)
            searchres = strfind(name,'c2c');
            if isempty(searchres)
                sampname = name;
            else
                etalname = name;
            end
        end
    end
    if ~isempty(strfind('Sao',sampleID))
        sampname = etalname;
    end
    if exist('sampname','var')==1
        sampdata = load([datadir,sampname]);
    else
        error('No sample data was loaded.')
    end
    if exist('etalname','var')==1
        etaldata = load([datadir,etalname]);
    else
        % if there's no c2c data Sthff are used
        for i = 1:size(filenamelist,1)
            name = filenamelist(i).name;
            searchres1 = strfind(name,['Sthff','_']);
            searchres2 = strfind(name,['Sthff','.']);
            searchres3 = strfind(name,['Sthff',' ']);
            searchres = [searchres1 searchres2 searchres3];
            clear searchres1 searchres2 searchres3;
            if ~isempty(searchres)
                searchres = strfind(name,'f2f');
                if ~isempty(searchres)
                    etalname = name;
                    fprintf(['WARNING. For ',sampleID,' no c2c data, Sthff c2c are used.\n']); 
                end
            end
        end
        etaldata = load([datadir,etalname]);
    end

    valuenum = max([size(sampdata.DataP,1),size(sampdata.DataS,1)]);
    tracenum = max([size(sampdata.DataP,2),size(sampdata.DataS,2)])+2;
    if strcmp(traceList(identif,1),sampleID);
        reprTraces = traceList(identif,2); reprTraces = reprTraces{1};
    else
        error('Wrong order for traceList. Please order like identifList.');
    end
	clear name searchres i;
    
    if ~exist('reprTracesP','var')==1
        reprTracesP = zeros(valuenum,size(identifList,1)*3);
        reprTracesS = zeros(valuenum,size(identifList,1)*3);
        reprTimeP = zeros(valuenum,size(identifList,1)*3);
        reprTimeS = zeros(valuenum,size(identifList,1)*3);
        reprMaxPos = zeros(size(identifList,1)*3,4);
        reprVpVs = zeros(size(identifList,1),2*3);
        reprQpQs = zeros(size(identifList,1),2*3);
        reprDataMIT = nan(9,size(identifList,1)*3);
    end

    %% Plotting information about pressure, temperature, height
    FigPTH = figure(figID.info+(identif-1)*length(fieldnames(figID)));
    set(FigPTH,'units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
    ax1 = subplot('Position',[0.05 0.66 0.90 0.30]); hold on;
    plot(3:tracenum,sampdata.DataMIT(5,:),'--b');
    [hAx,~,~] = plotyy(3:tracenum,sampdata.DataMIT(6,:),3:tracenum,mean([sampdata.DataMIT(3,:); sampdata.DataMIT(4,:)]));
    xlim(hAx(1),[3,tracenum]); xlim(hAx(2),[3,tracenum]);
    ylabel(hAx(1),'external pressure, atm'); ylabel(hAx(2),'pore pressure, ?');
    title(['Sample: ',sampleID]);
    legend('P_{axial}','P_{lateral}','P_{pore}');
    grid on;
    ax2 = subplot('Position',[0.05 0.36 0.90 0.25]); hold on;
    plot(3:tracenum,sampdata.DataMIT(2,:),'k');
    xlim([3,tracenum]);
    ylabel('temperature, {}^{o}C');
    legend('T_{inside camera}');
    grid on;
    ax3 = subplot('Position',[0.05 0.06 0.90 0.25]); hold on;
    plot(3:tracenum,sampdata.DataMIT(9,:),'k');
    xlim([3,tracenum]);
    xlabel('trace number'); ylabel('height, mm');
    legend('h');
    grid on;
    linkaxes([ax1,ax2,ax3],'x');
    clear ax1 ax2 ax3;
    if settSavePlot.saveFigsMAT
        savefig(FigPTH,[thisFolderName,'\outputs\PTH_',sampleID]);
    end
    if settSavePlot.saveFigsPNG
        set(FigPTH,'PaperPositionMode','auto');
        saveas(FigPTH,[thisFolderName,'\outputs\PTH_',sampleID],'png');
    end

    %% Load manual picks for current identif
    fileID = fopen(['data\fbs_hydrate_exp\manualPick_',sampleID,'_P.dat'],'r');
    if fileID < 0
        manualPickP = [];
    else
        manualPickP = fscanf(fileID,'%d %f',[2 Inf]);
        manualPickP = manualPickP';
        fclose(fileID);
    end
    fileID = fopen(['data\fbs_hydrate_exp\manualPick_',sampleID,'_S.dat'],'r');
    if fileID < 0
        manualPickS = [];
    else
        manualPickS = fscanf(fileID,'%d %f',[2 Inf]);
        manualPickS = manualPickS';
        fclose(fileID);
    end

    %% Prepare data and choose right etalon trace
    meanH   = mean(sampdata.DataMIT(9,:));
    % choose c2c trace
    heights = etaldata.DataMIT(9,:);
    pressures = etaldata.DataMIT(5,:);
    [~,posF2F] = max(pressures(heights > 3 & heights < 10)); % !!! changed from 3 to-1
    posF2F = posF2F + find(heights > 3 & heights < 10,1,'First') - 1; % !!! changed from 3 to-1
    if isempty(posF2F)
        fprintf('Please check etalon data. There is no heights > 3 mm and < 10 mm.\n'); 
        error('>>> STOP');
    end
    % choose aluminum sample trace
    heights = dataAl.DataMIT(10,:);
    [~,posAl] = sort(abs(heights - meanH));
    posAl = posAl(1,6); % choose sample with max pressure
    Time = [etaldata.Time(:,1) dataAl.Time(:,posAl) sampdata.Time(:,1)];
    clear meanPax meanH heights pressures;

    Velocities  = NaN(tracenum-2,2);
    Attenuation = NaN(tracenum-2,2);
    MaxPosInfo  = NaN(tracenum,4);
    AttenFromT  = NaN(tracenum-2,2);
    maxpos  = NaN(tracenum,3);
    quarPer = NaN(tracenum,2);
    for wavetype = 1:2
        switch wavetype
            case 1
                Data = [etaldata.DataP(:,posF2F) dataAl.DataP(:,posAl) sampdata.DataP(:,:)];
                manualPick = manualPickP;
            case 2
                Data = [etaldata.DataS(:,posF2F) dataAl.DataS(:,posAl) sampdata.DataS(:,:)];
                manualPick = manualPickS;
        end
        badIDs = identifList(identif,wavetype+1);
        badIDs = badIDs{1};
        
        %% Change NaNs to 0
        Data(isnan(Data)) = 0;
        %% Correcting shift up
        Data = Data - repmat(mean(Data),size(Data,1),1);
        %% Noise filtering
        if settV.noiseFilter
            for trace = 1:tracenum
                % choose time (f2f, aluminum or sample)
                if trace > 2
                    timeID = 3;
                else
                    timeID = trace;
                end
                df = 1/(Time(end,timeID)-Time(1,timeID));
                freq = 0:df:df*(size(Time,1)-1);
                filt = ones(1,size(freq,2));
                pos1 = floor(settV.filtFreqRange(1,1)/df)+1;
                pos2 = ceil(settV.filtFreqRange(1,2)/df)+1;
                pos3 = floor(settV.filtFreqRange(1,3)/df)+1;
                pos4 = ceil(settV.filtFreqRange(1,4)/df)+1;
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

        %% Find first max position using the windows
        for trace = 1:tracenum
            % choose time (f2f, aluminum or sample)
            if trace > 2
                timeID = 3;
            else
                timeID = trace;
            end
            % skip bad traces, IDs from identifList
            if any(badIDs == trace)
                maxpos(trace,:) = [NaN NaN NaN];
                continue;
            end
            % accepting manual picks
            temptime = NaN;
            if ~isnan(manualPick)
                pos = find(manualPick(:,1) == trace);
                if ~isempty(pos)
                    temptime = manualPick(pos,2);
                    maxpos(trace,3) = 1;
                end
            end
            % if there is no manual pick
            if isnan(temptime)
                % find max positions for polynomial interpolation
                findval = settV.firstVal*max(Data(:,trace));
                rangepos = find(Data(:,trace)>findval,1,'first');
                if isempty(rangepos)
                    maxpos(trace,:) = [NaN NaN NaN];
                    continue;
                end
                % cut signal using maxWind and find max
                cutSig = cutSignal([Time(:,timeID) Data(:,trace)], ...
                                   [Time(rangepos,timeID),Time(rangepos,timeID)+settV.maxWind]);
                [val,pos] = max(cutSig(:,2));
                temptime = Time(rangepos+pos-1,timeID);
                maxpos(trace,3) = 0;
            end
            % cut signal using intWind and interpolate
            cutSig = cutSignal([Time(:,timeID) Data(:,trace)], ...
                               [temptime-settV.intWind/2,temptime+settV.intWind/2]);
            inttime = cutSig(1,1):settV.polyfitdisc:cutSig(end,1);
            ws = warning('off','all');
            intpoly = polyfit(cutSig(:,1),cutSig(:,2),settV.polyfitdegr);
            warning(ws);
            [maxpos(trace,2),pos] = max(polyval(intpoly,inttime));
            maxpos(trace,1) = inttime(1,pos);
        end
        clear timeID badIDs findval val rangepos pos temptime inttime intpoly trace ws cutSig;

        %% Velocity calculation
        smplsize = (1.0e-3*sampdata.DataMIT(9,:)-7.12e-3)'; % sample heights, m
        trvltime = maxpos(:,1)-maxpos(1,1); trvltime(1:2) = []; % traveltimes, s
        % different for some THF experiments
        if ~isempty(strfind('Sthfa_Sthfb_Sthfc_Sthfd_Sthfe_Sthff_Sthfg_Sthfh_Sthfi',sampleID))
            heights = etaldata.DataMIT(9,:);
            pos = find(heights > 3 & heights < 10,1,'Last');
            if wavetype==1
                trvltime = maxpos(:,1)-( 6.672e-6 - 1.0e-3*(7.16-heights(1,pos))/5950 );
            else
                trvltime = maxpos(:,1)-( 11.76e-6 - 1.0e-3*(7.16-heights(1,pos))/3230 );
            end
            trvltime(1:2) = [];
            clear heights pos;
        end
        Velocities(:,wavetype) = smplsize./trvltime; % P, S velocities, m/s
        clear smplsize trvltime;
        if wavetype==1
            MaxPosInfo(:,1:2) = maxpos(:,1:2);
        else
            MaxPosInfo(:,3:4) = maxpos(:,1:2);
        end

        %% Plotting traces
        if settSavePlot.plotTraces
            if wavetype==1
                ID = figID.tracesP+(identif-1)*length(fieldnames(figID));
                pickspath = ['data\fbs_hydrate_exp\manualPick_',sampleID,'_P.dat'];
            else
                ID = figID.tracesS+(identif-1)*length(fieldnames(figID));
                pickspath = ['data\fbs_hydrate_exp\manualPick_',sampleID,'_S.dat'];
            end
            FigTrace = figure(ID);
            clear ID;
            set(FigTrace,'units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
            subplot('Position',[0.05 0.06 0.90 0.88]); hold on;
            if settSavePlot.plotTracesNorm
                maxval = max(abs(Data));
            else
                maxval = [max(abs(Data(:,1))) max(abs(Data(:,2))) ...
                          repmat(max(max(abs(Data(:,3:end)))),1,size(Data,2)-2)];
            end
            tempData = repmat(1:size(Data,2),size(Data,1),1) + Data./repmat(maxval,size(Data,1),1)/2;
            plot([Time(:,1) Time(:,2) repmat(Time(:,3),1,size(Data,2)-2)],tempData,'k');
            tempData = (1:size(maxpos,1))' + maxpos(:,2)./maxval'/2;
            plot(maxpos(:,1),tempData,'+r');
            pos = find(maxpos(:,3));
            plot(maxpos(pos,1),tempData(pos,1),'+g');
            clear maxval tempData pos;
            xlim([0,max([Time(end,1),Time(end,2),Time(end,3)])]);
            ylim([0,20]);
            LinePlotExplorer_polyfit(pickspath);
            % figure settings for traces
            set(gca,'Ytick',1:tracenum); grid on; xlabel('times, s'); ylabel('trace number');
            clear pickspath;
        end
        
        %% Saving results in file
        switch wavetype
            case 1
                title(['Sample: ',sampleID,',  P waves,  ',num2str(tracenum),' traces']);
                dlmwrite([thisFolderName,'\outputs\',sampleID,'_P_fbs.dat'],maxpos,'\t');
                if settSavePlot.saveFigsMAT
                    savefig(FigTrace,[thisFolderName,'\outputs\P_traces_',sampleID]);
                end
                if settSavePlot.saveFigsPNG && settSavePlot.plotTraces
                    set(FigTrace,'PaperPositionMode','auto');
                    saveas(FigTrace,[thisFolderName,'\outputs\P_traces_',sampleID],'png');
                end
            case 2
                title(['Sample: ',sampleID,',  S waves,  ',num2str(tracenum),' traces']);
                dlmwrite([thisFolderName,'\outputs\',sampleID,'_S_fbs.dat'],maxpos,'\t');
                if settSavePlot.saveFigsMAT
                    savefig(FigTrace,[thisFolderName,'\outputs\S_traces_',sampleID]);
                end
                if settSavePlot.saveFigsPNG && settSavePlot.plotTraces
                    set(FigTrace,'PaperPositionMode','auto');
                    saveas(FigTrace,[thisFolderName,'\outputs\S_traces_',sampleID],'png');
                end
        end
        
        %% Q estimation
        if settQ.Qestim
            switch wavetype
                case 1
                    numOf0Cross = settQ.numOf0CrossP;
                case 2
                    numOf0Cross = settQ.numOf0CrossS;
            end
            cutSign = NaN(valuenum,tracenum);
            cutTime = NaN(valuenum,tracenum);
            % skip first (f2f) trace
            for trace = 2:tracenum
                % choose time (f2f, aluminum or sample)
                if trace > 2
                    timeID = 3;
                else
                    timeID = trace;
                end
                % skip bad traces, IDs from identifList
                badIDs = identifList(identif,wavetype+1); badIDs = badIDs{1};
                if any(badIDs == trace)
                    continue;
                end
                % search for zero cross
                time = Time(:,timeID);
                data = Data(:,trace);
                data = data(time>maxpos(trace,1));
                time = time(time>maxpos(trace,1));
                [~,t0] = crossing(data,time,0,'linear');
                quarPer(trace,wavetype) = t0(1,1)-maxpos(trace,1);
                % cutting signal
                windstr = {'timewindP', 'timewindS'};
                if isfield(settQ,windstr{wavetype})
                    if wavetype==1
                        timewind = settQ.timewindP;
                    else
                        timewind = settQ.timewindS;
                    end
                    [~,pos] = min(abs(t0 - ...
                              (maxpos(trace,1)-timewind(1,1)+timewind(1,2)) ));
                    temp = cutSignal([Time(:,timeID) Data(:,trace)], ...
                        [maxpos(trace,1)-timewind(1,1),t0(1,pos)]);
                    cutTime(1:size(temp,1),trace) = temp(:,1);
                    cutSign(1:size(temp,1),trace) = temp(:,2);
                    clear pos;
                else
                    temp = cutSignal([Time(:,timeID) Data(:,trace)], ...
                        [maxpos(trace,1)-2*quarPer(trace,wavetype),t0(1,numOf0Cross)]);
                    cutTime(1:size(temp,1),trace) = temp(:,1);
                    cutSign(1:size(temp,1),trace) = temp(:,2);
                end
            end
            clear timeID t0 time data temp;
            % skip NaNs (aluminum)
            alumInd = 2;
            pos = find(isnan(cutTime(:,alumInd)),1,'first');
            if pos == 1
                error('No aluminum signal after cutting.');
            end
            sigCutBef = [cutTime(1:pos-1,alumInd),cutSign(1:pos-1,alumInd)];
            % adding window smoothing (aluminum)
            sigCutBef = windows(sigCutBef, settQ.dtw, 'blackman');
            clear alumInd;

            % preparing figures
            FigSigAndSp = figure(figID.SigAndSp+(identif-1)*length(fieldnames(figID)));
            set(FigSigAndSp,'units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
            FigLnSp = figure(figID.lnSp+(identif-1)*length(fieldnames(figID)));
            set(FigLnSp,'units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
            % skipped first two traces (f2f and Al)
            for trace = 3:tracenum
                % skip bad traces, IDs from identifList
                badIDs = identifList(identif,wavetype+1); badIDs = badIDs{1};
                if any(badIDs == trace)
                    Attenuation(trace-2,wavetype) = NaN;
                    continue;
                end
                % skip NaNs
                pos = find(isnan(cutTime(:,trace)),1,'first');
                if pos == 1
                    Attenuation(trace-2,wavetype) = NaN;
                    continue;
                end
                sigCutAft = [cutTime(1:pos-1,trace),cutSign(1:pos-1,trace)];
                % adding window smoothing
                sigCutAft = windows(sigCutAft, settQ.dtw, 'blackman');
                % adding zeros for equal time length
                maxt = max([sigCutBef(end,1)-sigCutBef(1,1) sigCutAft(end,1)-sigCutAft(1,1)]);
                sigCutBef0 = addZeros(sigCutBef, maxt);
                sigCutAft0 = addZeros(sigCutAft, maxt);
                % calculating amplitude spectrums
                spAmpBef = spectrumAmp(sigCutBef0);
                spAmpAft = spectrumAmp(sigCutAft0);
                [spAmpBef, spAmpAft] = checkSpectrPair(spAmpBef, spAmpAft);
                % plotting
                if settSavePlot.plotSigAndSp
                    if settSavePlot.plotReprOnly
                        if any(trace == reprTraces)
                            figure(figID.SigAndSp+(identif-1)*length(fieldnames(figID)));
                            subplot('Position',[(wavetype-1)*0.5+0.05 0.55 0.43 0.41]); hold on;
                            plot(sigCutAft(:,1)-sigCutAft(1,1),sigCutAft(:,2));
                            subplot('Position',[(wavetype-1)*0.5+0.05 0.06 0.43 0.41]); hold on;
                            plot(spAmpAft(:,1),spAmpAft(:,2));
                        end
                    else
                        if any(trace == 2:settQ.plotSkipNum+1:tracenum-1)
                            figure(figID.SigAndSp+(identif-1)*length(fieldnames(figID)));
                            subplot('Position',[(wavetype-1)*0.5+0.05 0.55 0.43 0.41]); hold on;
                            plot(sigCutAft(:,1)-sigCutAft(1,1),sigCutAft(:,2));
                            subplot('Position',[(wavetype-1)*0.5+0.05 0.06 0.43 0.41]); hold on;
                            plot(spAmpAft(:,1),spAmpAft(:,2));
                        end
                    end
                end

                switch wavetype
                    case 1
                        freqRange = settQ.freqRangeP;
                    case 2
                        freqRange = settQ.freqRangeS;
                end
                % calculating Q
                [invQ,invQerr,LnSpOut,linfit] = ...
                    QfromSpRatio(spAmpBef,spAmpAft,freqRange,maxpos(trace,1)-maxpos(1,1),'off','','','');
                if invQ < -1 || invQ > 1
                    invQ = NaN;
                end
                Attenuation(trace-2,wavetype) = invQ;
                % plotting
                if settSavePlot.plotLnSp
                    if settSavePlot.plotReprOnly
                        if any(trace == reprTraces)
                            figure(figID.lnSp+(identif-1)*length(fieldnames(figID)));
                            subplot('Position',[0.05 (1-wavetype)*0.5+0.56 0.93 0.41]); hold on;
                            plot(freqRange,polyval(linfit,freqRange),'k','LineWidth',1.0);
                            plot(LnSpOut(:,1),LnSpOut(:,2),'-*');
                        end
                    else
                        if any(trace == 3:settQ.plotSkipNum+1:tracenum)
                            figure(figID.lnSp+(identif-1)*length(fieldnames(figID)));
                            subplot('Position',[0.05 (1-wavetype)*0.5+0.56 0.93 0.41]); hold on;
                            plot(freqRange,polyval(linfit,freqRange),'k','LineWidth',1.0);
                            plot(LnSpOut(:,1),LnSpOut(:,2),'-*');
                        end
                    end
                end
                
                % Estimation attenuation from quarter of period
                AttenFromT(trace-2,wavetype) = 4*(quarPer(trace,wavetype)-quarPer(2,wavetype))/ ...
                                               (maxpos(trace,1)-maxpos(1,1));
            end
            % add signal and spectrum for aluminum
            if settSavePlot.plotSigAndSp
                figure(figID.SigAndSp+(identif-1)*length(fieldnames(figID)));
                subplot('Position',[(wavetype-1)*0.5+0.05 0.55 0.43 0.41]); hold on;
                temp = plot(sigCutBef(:,1)-sigCutBef(1,1),sigCutBef(:,2),'k','LineWidth',2);
                legend(temp,'aluminum');
                subplot('Position',[(wavetype-1)*0.5+0.05 0.06 0.43 0.41]); hold on;
                temp = plot(spAmpBef(:,1),spAmpBef(:,2),'k','LineWidth',2);
                legend(temp,'aluminum');
                clear temp;
                % figure settings for signals and spectrums
                subplot('Position',[0.05 0.55 0.43 0.41]); title(['Sample: ',sampleID,',  P waves, signals']);
                grid on; xlabel('time, s');
                subplot('Position',[0.55 0.55 0.43 0.41]); title(['Sample: ',sampleID,',  S waves, signals']);
                grid on; xlabel('time, s');
                subplot('Position',[0.05 0.06 0.43 0.41]); title(['Sample: ',sampleID,',  P waves, spectrums']);
                xlim([0,settQ.plotFreqMax]); grid on; xlabel('frequency, Hz'); ylabel('amplitude spectrum');
                subplot('Position',[0.55 0.06 0.43 0.41]); title(['Sample: ',sampleID,',  S waves, spectrums']);
                xlim([0,settQ.plotFreqMax]); grid on; xlabel('frequency, Hz'); ylabel('amplitude spectrum');
                if settSavePlot.saveFigsMAT
                    savefig(figure(figID.SigAndSp+(identif-1)*length(fieldnames(figID))),[thisFolderName,'\outputs\SigAndSp_',sampleID]);
                end
                if settSavePlot.saveFigsPNG && settSavePlot.plotSigAndSp
                    set(figure(figID.SigAndSp+(identif-1)*length(fieldnames(figID))),'PaperPositionMode','auto');
                    saveas(figure(figID.SigAndSp+(identif-1)*length(fieldnames(figID))),[thisFolderName,'\outputs\SigAndSp_',sampleID],'png');
                end
            end
            if settSavePlot.plotLnSp
                figure(figID.lnSp+(identif-1)*length(fieldnames(figID)));
                subplot('Position',[0.05 0.56 0.93 0.41]); title(['Sample: ',sampleID,',  P waves, Ln(Sp)']);
                xlim([0,settQ.plotFreqMax]); grid on; xlabel('frequency, Hz');
                subplot('Position',[0.05 0.06 0.93 0.41]); title(['Sample: ',sampleID,',  S waves, Ln(Sp)']);
                xlim([0,settQ.plotFreqMax]); grid on; xlabel('frequency, Hz');
                if settSavePlot.saveFigsMAT
                    savefig(figure(figID.lnSp+(identif-1)*length(fieldnames(figID))),[thisFolderName,'\outputs\LnSp_',sampleID]);
                end
                if settSavePlot.saveFigsPNG && settSavePlot.plotLnSp
                    set(figure(figID.lnSp+(identif-1)*length(fieldnames(figID))),'PaperPositionMode','auto');
                    saveas(figure(figID.lnSp+(identif-1)*length(fieldnames(figID))),[thisFolderName,'\outputs\LnSp_',sampleID],'png');
                end
            end

            clear trace pos cutTime cutSign sigCutBef sigCutAft sigCutBef0 sigCutAft0;
            clear spAmpBef spAmpAft freqRange invQ invQerr LnSpOut linfit maxt numOf0Cross;
        end
        
        %% Collecting representative traces
        if wavetype==1
            reprTimeP(:,identif)                       = Time(:,3);
            reprTimeP(:,identif+size(identifList,1))   = Time(:,3);
            reprTimeP(:,identif+size(identifList,1)*2) = Time(:,3);
            if reprTraces(1,1)~= 0
                reprTracesP(:,identif) = Data(:,reprTraces(1,1));
            end
            if reprTraces(1,2)~= 0
                reprTracesP(:,identif+size(identifList,1)) = Data(:,reprTraces(1,2));
            end
            if reprTraces(1,3)~= 0
                reprTracesP(:,identif+size(identifList,1)*2) = Data(:,reprTraces(1,3));
            end
        else
            reprTimeS(:,identif)                       = Time(:,3);
            reprTimeS(:,identif+size(identifList,1))   = Time(:,3);
            reprTimeS(:,identif+size(identifList,1)*2) = Time(:,3);
            if reprTraces(1,1)~= 0
                reprTracesS(:,identif) = Data(:,reprTraces(1,1));
                reprMaxPos(identif,:) = MaxPosInfo(reprTraces(1,1),:);
                reprVpVs(identif,[1 4]) = Velocities(reprTraces(1,1)-2,:);
                reprQpQs(identif,[1 4]) = Attenuation(reprTraces(1,1)-2,:);
                reprDataMIT(:,identif) = sampdata.DataMIT(1:9,reprTraces(1,1)-2);
            end
            if reprTraces(1,2)~= 0
                reprTracesS(:,identif+size(identifList,1)) = Data(:,reprTraces(1,2));
                reprMaxPos(identif+size(identifList,1),:) = MaxPosInfo(reprTraces(1,2),:);
                reprVpVs(identif,[2 5]) = Velocities(reprTraces(1,2)-2,:);
                reprQpQs(identif,[2 5]) = Attenuation(reprTraces(1,2)-2,:);
                reprDataMIT(:,identif+size(identifList,1)) = sampdata.DataMIT(1:9,reprTraces(1,2)-2);
            end
            if reprTraces(1,3)~= 0
                reprTracesS(:,identif+size(identifList,1)*2) = Data(:,reprTraces(1,3));
                reprMaxPos(identif+size(identifList,1)*2,:) = MaxPosInfo(reprTraces(1,3),:);
                reprVpVs(identif,[3 6]) = Velocities(reprTraces(1,3)-2,:);
                reprQpQs(identif,[3 6]) = Attenuation(reprTraces(1,3)-2,:);
                reprDataMIT(:,identif+size(identifList,1)*2) = sampdata.DataMIT(1:9,reprTraces(1,3)-2);
            end
        end
    end
    clear posF2F posAl;

    %% Saving results in file
    dlmwrite([thisFolderName,'\outputs\',sampleID,'_VpVs.dat'],Velocities,'\t');
    dlmwrite([thisFolderName,'\outputs\',sampleID,'_QpQs.dat'],Attenuation,'\t');
    dlmwrite([thisFolderName,'\outputs\',sampleID,'_quarPer.dat'],quarPer,'\t');

    %% Plotting results (velocities, attenuation, pressure, temperature)
    FigVQ = figure(figID.results+(identif-1)*length(fieldnames(figID)));
    set(FigVQ,'units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
    ax1 = subplot('Position',[0.05 0.66 0.90 0.30]); hold on;
	plot(3:tracenum,Velocities(:,1),'b');
    plot(3:tracenum,Velocities(:,2),'r');
    xlim([3,tracenum]);
    ylabel('velocity, m/s');
    title(['Sample: ',sampleID]); 
    legend('V_P','V_S');
    grid on;
    ax2 = subplot('Position',[0.05 0.36 0.90 0.25]); hold on;
    plot(3:tracenum,Attenuation(:,1),'b');
    plot(3:tracenum,Attenuation(:,2),'r');
    xlim([3,tracenum]);
    ylabel('attenuation');
    legend('Q_P^{-1}','Q_S^{-1}');
    grid on;
    ax3 = subplot('Position',[0.05 0.06 0.90 0.25]);
    [hAx,~,~] = plotyy(3:tracenum,sampdata.DataMIT(2,:),3:tracenum,sampdata.DataMIT(5,:));
    hold(hAx(2)); plot(hAx(2),3:tracenum,sampdata.DataMIT(6,:));
    temp = [sampdata.DataMIT(5,:) sampdata.DataMIT(6,:)];
    temp = [floor(min(temp)-0.1*(max(temp)-min(temp))),ceil(max(temp)+0.1*(max(temp)-min(temp)))];
    temp(1,3) = ceil((temp(1,2)-temp(1,1))/(numel(get(hAx(1),'YTick'))-1));
    temp(1,4) = temp(1,1)+temp(1,3)*(numel(get(hAx(1),'YTick'))-1);
    set(hAx(2),'YTick',temp(1,1):temp(1,3):temp(1,4));
    ylim(hAx(2),[temp(1,1),temp(1,4)]);
    xlim(hAx(1),[3,tracenum]); xlim(hAx(2),[3,tracenum]);
    ylabel(hAx(1),'temperature, {}^{o}C'); ylabel(hAx(2),'pressure, atm');
    legend('T_{inside}','P_{axial}','P_{lateral}');
    grid on;
    linkaxes([ax1,ax2,ax3,hAx(2)],'x');
    clear ax1 ax2 ax3;
    if settSavePlot.saveFigsMAT
        savefig(FigVQ,[thisFolderName,'\outputs\VpVsQpQs_',sampleID]);
    end
    if settSavePlot.saveFigsPNG
        set(FigVQ,'PaperPositionMode','auto');
        saveas(FigVQ,[thisFolderName,'\outputs\VpVsQpQs_',sampleID],'png');
    end
    
    %% Plotting results velocities, attenuation dependence from temperature
    FigVQT = figure(figID.resVQfT+(identif-1)*length(fieldnames(figID)));
    set(FigVQT,'units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
    ax1 = subplot('Position',[0.06 0.06 0.43 0.90]); hold on;
	plot(sampdata.DataMIT(2,:),Velocities(:,1),'b');
    plot(sampdata.DataMIT(2,:),Velocities(:,2),'r');
    xlabel('temperature, C^o');
    ylabel('velocity, m/s');
    legend('V_P','V_S');
    grid on; title(['Sample: ',sampleID]);
    ax2 = subplot('Position',[0.55 0.06 0.43 0.90]); hold on;
    plot(sampdata.DataMIT(2,:),Attenuation(:,1),'b');
    plot(sampdata.DataMIT(2,:),Attenuation(:,2),'r');
    xlabel('temperature, C^o');
    ylabel('attenuation');
    legend('Q_P^{-1}','Q_S^{-1}');
    grid on; title(['Sample: ',sampleID]);
    linkaxes([ax1,ax2],'x');
    clear ax1 ax2 ax3;
    if settSavePlot.saveFigsMAT
        savefig(FigVQT,[thisFolderName,'\outputs\VQfromT_',sampleID]);
    end
    if settSavePlot.saveFigsPNG
        set(FigVQT,'PaperPositionMode','auto');
        saveas(FigVQT,[thisFolderName,'\outputs\VQfromT_',sampleID],'png');
    end
    
    %% Saving results in xls file
    if settSavePlot.saveXLS
        fprintf('Saving results in XLS file, please wait.\n');
        path = [thisFolderName,'\outputs\',sampname,'.xls'];
        if exist(path,'file')==2
            delete(path);
        end
        ws = warning('off','all');
        xlsdata = [{'trace number'} {'time from start, hours'} {'Vp, m/s'} {'Vs, m/s'} {'Qp'} {'Qs'}; ...
                   num2cell((1:tracenum-2)') num2cell((sampdata.DataMIT(1,:)/3600)') ...
                   num2cell(Velocities) num2cell(Attenuation)];
        xlswrite(path,xlsdata,'results');
        if size(sampdata.DataMIT,1)==16
            xlsdata = [{'trace number'} {'time from start, s'} {'temperature inside camera, C'} ...
                       {'pore pressurea at the inlet, mV'} {'pore pressure at the outlet, mV'} ...
                       {'axial pressure, atm'} {'radial pressure, atm'} ...
                       {'-not used-'} {'room temperature, C'} {'sample height, mm'} ...
                       {'programmed temperature, C'} {'process temperature, C'} ...
                       {'temperature inside thermostat, C'} {''} {''} {''} {''}; ...
                       num2cell((1:tracenum-2)') num2cell((sampdata.DataMIT)')];
        else
            if size(sampdata.DataMIT,1)==12
                xlsdata = [{'trace number'} {'time from start, s'} {'temperature inside camera, C'} ...
                           {'pore pressurea at the inlet, mV'} {'pore pressure at the outlet, mV'} ...
                           {'axial pressure, atm'} {'radial pressure, atm'} ...
                           {'-not used-'} {'room temperature, C'} {'sample height, mm'} ...
                           {'programmed temperature, C'} {'process temperature, C'} ...
                           {'temperature inside thermostat, C'}; ...
                           num2cell((1:tracenum-2)') num2cell((sampdata.DataMIT)')];
            else
                xlsdata = [{'trace number'} {'time from start, s'} {'temperature inside camera, C'} ...
                           {'pore pressurea at the inlet, mV'} {'pore pressure at the outlet, mV'} ...
                           {'axial pressure, atm'} {'radial pressure, atm'} ...
                           {'-not used-'} {'room temperature, C'} {'sample height, mm'}; ...
                           num2cell((1:tracenum-2)') num2cell((sampdata.DataMIT)')];
            end
        end
        xlswrite(path,xlsdata,'DataMIT');
        xlsdata = [{'trace number'} {'max positions for P, s'} {'maxpos amplitudes for P'} ...
                   {'max positions for S, s'} {'maxpos amplitudes for S'} ...
                   {'quarter of period for P, s'} {'quarter of period for S, s'} ...
                   {'attenuation from T/4 for P'} {'attenuation from T/4 for S'}; ...
                   [{'f2f'; 'Al'}; num2cell((1:tracenum-2)')] num2cell(MaxPosInfo) ...
                   num2cell(quarPer) num2cell([NaN(2,2); AttenFromT])];
        xlswrite(path,xlsdata,'fb_amp_quarT_Q');
        xlsdata = [{'file with etalon data:'} {etalname}; {'file with sample data:'} {sampname}];
        xlswrite(path,xlsdata,'used files');
        warning(ws);

        %% Plotting in excel file
        if settSavePlot.plotXLS
            fprintf('Plotting in XLS file.\n');
            if exist(path,'file')==2
                e = actxserver('excel.application');
                eWs = e.Workbooks;
                eW = eWs.Open(path);
                %e.Visible = 1;
                eS  = eW.Sheets.Item('results');
                eS2 = eW.Sheets.Item('DataMIT');
                eS3 = eW.Sheets.Item('fb_amp_quarT_Q');
                % plot velocity data
                eCO1 = eS.ChartObjects.Add(300, 10, 600, 300);
                eC1 = eCO1.Chart;
                eC1.SeriesCollection.NewSeries;
                eC1.SeriesCollection(1).Value = eS.Range(['C2:C',num2str(tracenum-1)]);
                eC1.SeriesCollection(1).XValue = eS.Range(['B2:B',num2str(tracenum-1)]);
                eC1.SeriesCollection(1).Name = 'Vp';
                eC1.SeriesCollection(1).Format.Line.ForeColor.RGB = 255*65536 + 0*256 + 0;
                eC1.SeriesCollection.NewSeries;
                eC1.SeriesCollection(2).Value = eS.Range(['D2:D',num2str(tracenum-1)]);
                eC1.SeriesCollection(2).XValue = eS.Range(['B2:B',num2str(tracenum-1)]);
                eC1.SeriesCollection(2).Name = 'Vs';
                eC1.SeriesCollection(2).Format.Line.ForeColor.RGB = 0*65536 + 0*256 + 255;
                eCO1.Chart.ChartType = 'xlXYScatterLinesNoMarkers';
                eCO1.Chart.Legend.Position = 'xlLegendPositionBottom';
                eCO1.Chart.HasTitle = true;
                eCO1.Chart.ChartTitle.Text = 'Velocities';
                eCO1.Chart.Axes(1).HasTitle = true;
                eCO1.Chart.Axes(1).AxisTitle.Caption = 'time from start, hours';
                eCO1.Chart.Axes(1).AxisTitle.Characters.Font.Bold = false;
                eCO1.Chart.Axes(1).HasMajorGridlines = true;
                % plot temperature data
                eCO2 = eS.ChartObjects.Add(300, 320, 600, 300);
                eC2 = eCO2.Chart;
                eC2.SeriesCollection.NewSeries;
                eC2.SeriesCollection(1).Value = eS2.Range(['C2:C',num2str(tracenum-1)]);
                eC2.SeriesCollection(1).XValue = eS.Range(['B2:B',num2str(tracenum-1)]);
                eC2.SeriesCollection(1).Name = 'T';
                eC2.SeriesCollection(1).Format.Line.ForeColor.RGB = 0*65536 + 255*256 + 0;
                eCO2.Chart.ChartType = 'xlXYScatterLinesNoMarkers';
                eCO2.Chart.Legend.Position = 'xlLegendPositionBottom';
                eCO2.Chart.HasTitle = true;
                eCO2.Chart.ChartTitle.Text = 'Temperature';
                eCO2.Chart.Axes(1).HasTitle = true;
                eCO2.Chart.Axes(1).AxisTitle.Caption = 'time from start, hours';
                eCO2.Chart.Axes(1).AxisTitle.Characters.Font.Bold = false;
                eCO2.Chart.Axes(1).HasMajorGridlines = true;
                % plot attenuation data
                eCO3 = eS.ChartObjects.Add(300, 630, 600, 300);
                eC3 = eCO3.Chart;
                eC3.SeriesCollection.NewSeries;
                eC3.SeriesCollection(1).Value = eS.Range(['E2:E',num2str(tracenum-1)]);
                eC3.SeriesCollection(1).XValue = eS.Range(['B2:B',num2str(tracenum-1)]);
                eC3.SeriesCollection(1).Name = 'Qp';
                eC3.SeriesCollection(1).Format.Line.ForeColor.RGB = 255*65536 + 0*256 + 0;
                eC3.SeriesCollection.NewSeries;
                eC3.SeriesCollection(2).Value = eS.Range(['F2:F',num2str(tracenum-1)]);
                eC3.SeriesCollection(2).XValue = eS.Range(['B2:B',num2str(tracenum-1)]);
                eC3.SeriesCollection(2).Name = 'Qs';
                eC3.SeriesCollection(2).Format.Line.ForeColor.RGB = 0*65536 + 0*256 + 255;
                eCO3.Chart.ChartType = 'xlXYScatterLinesNoMarkers';
                eCO3.Chart.Legend.Position = 'xlLegendPositionBottom';
                eCO3.Chart.HasTitle = true;
                eCO3.Chart.ChartTitle.Text = 'Attenuation';
                eCO3.Chart.Axes(1).HasTitle = true;
                eCO3.Chart.Axes(1).AxisTitle.Caption = 'time from start, hours';
                eCO3.Chart.Axes(1).AxisTitle.Characters.Font.Bold = false;
                eCO3.Chart.Axes(1).HasMajorGridlines = true;
                % plot additional attenuation data
                eCO4 = eS3.ChartObjects.Add(450, 10, 800, 500);
                eC4 = eCO4.Chart;
                eC4.SeriesCollection.NewSeries;
                eC4.SeriesCollection(1).Value = eS.Range(['E2:E',num2str(tracenum-1)]);
                eC4.SeriesCollection(1).XValue = eS.Range(['B2:B',num2str(tracenum-1)]);
                eC4.SeriesCollection(1).Name = 'Qp';
                eC4.SeriesCollection(1).Format.Line.ForeColor.RGB = 255*65536 + 0*256 + 0;
                eC4.SeriesCollection.NewSeries;
                eC4.SeriesCollection(2).Value = eS.Range(['F2:F',num2str(tracenum-1)]);
                eC4.SeriesCollection(2).XValue = eS.Range(['B2:B',num2str(tracenum-1)]);
                eC4.SeriesCollection(2).Name = 'Qs';
                eC4.SeriesCollection(2).Format.Line.ForeColor.RGB = 0*65536 + 0*256 + 255;
                eC4.SeriesCollection.NewSeries;
                eC4.SeriesCollection(3).Value = eS3.Range(['H4:H',num2str(tracenum-1)]);
                eC4.SeriesCollection(3).XValue = eS.Range(['B2:B',num2str(tracenum-1)]);
                eC4.SeriesCollection(3).Name = 'Qp from T/4';
                eC4.SeriesCollection(3).Format.Line.ForeColor.RGB = 255*65536 + 255*256 + 0;
                eC4.SeriesCollection.NewSeries;
                eC4.SeriesCollection(4).Value = eS3.Range(['I4:I',num2str(tracenum-1)]);
                eC4.SeriesCollection(4).XValue = eS.Range(['B2:B',num2str(tracenum-1)]);
                eC4.SeriesCollection(4).Name = 'Qs from T/4';
                eC4.SeriesCollection(4).Format.Line.ForeColor.RGB = 0*65536 + 165*256 + 255;
                eCO4.Chart.ChartType = 'xlXYScatterLinesNoMarkers';
                eCO4.Chart.Legend.Position = 'xlLegendPositionBottom';
                eCO4.Chart.HasTitle = true;
                eCO4.Chart.ChartTitle.Text = 'Attenuation';
                eCO4.Chart.Axes(1).HasTitle = true;
                eCO4.Chart.Axes(1).AxisTitle.Caption = 'time from start, hours';
                eCO4.Chart.Axes(1).AxisTitle.Characters.Font.Bold = false;
                eCO4.Chart.Axes(1).HasMajorGridlines = true;
                % save and close
                eW.Save;
                clear eCO1 eCO2 eC1 eC2 eS;
                eW.Close; eWs.Close; clear eW eWs;
                e.Quit; clear e;
                if size(identifList,1)>1 && identif<size(identifList,1)
                    pause(10);
                end
                if settSavePlot.saveSmpForXLS
                    set(FigVQ,'PaperPositionMode','auto');
                    saveas(FigVQ,[path(1:end-4),'.png'],'png');
                end
            else
                error('\nNo file:\n%s\nCan not plot in xls file.\n',path);
            end
        end
        clear path xlsdata;
    end
    if size(identifList,1) > 3
        close all;
    end
end

%% Plotting represantative traces
% for wavetype = 1:2
%     FigReprTrace = figure;
%     set(FigReprTrace,'units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
%     subplot('Position',[0.05 0.06 0.90 0.88]); hold on;
%     if wavetype==1
%         Time = reprTimeP;
%         Data = reprTracesP;
%     else
%         Time = reprTimeS;
%         Data = reprTracesS;
%     end
%     maxval = [max(abs(Data(:,1))) max(abs(Data(:,2))) ...
%               repmat(max(max(abs(Data(:,3:end)))),1,size(Data,2)-2)];
%     tempData = repmat(1:size(Data,2),size(Data,1),1) + Data./repmat(maxval,size(Data,1),1)/2;
%     plot(Time,tempData,'k');
%     tempData = (1:size(reprMaxPos,1))' + reprMaxPos(:,wavetype*2)./maxval'/2;
%     plot(reprMaxPos(:,wavetype*2-1),tempData,'+r');
%     clear maxval tempData pos;
%     xlim([0,max(Time(end,:))]);
%     ylim([0,20]);
%     LinePlotExplorer();
%     % figure settings for traces
%     set(gca,'Ytick',1:size(reprMaxPos,1)); set(gca,'YtickLabel',identifList(:,1));
%     grid on; xlabel('times, s'); ylabel('trace number');
%     if wavetype==1
%         title('P waves, real amplitudes');
%     else
%         title('S waves, real amplitudes');
%     end
% end
% for wavetype = 1:2
%     FigReprTrace = figure;
%     set(FigReprTrace,'units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
%     subplot('Position',[0.05 0.06 0.90 0.88]); hold on;
%     if wavetype==1
%         Time = reprTimeP;
%         Data = reprTracesP;
%     else
%         Time = reprTimeS;
%         Data = reprTracesS;
%     end
%     maxval = max(abs(Data));
%     tempData = repmat(1:size(Data,2),size(Data,1),1) + Data./repmat(maxval,size(Data,1),1)/2;
%     plot(Time,tempData,'k');
%     tempData = (1:size(reprMaxPos,1))' + reprMaxPos(:,wavetype*2)./maxval'/2;
%     plot(reprMaxPos(:,wavetype*2-1),tempData,'+r');
%     clear maxval tempData pos;
%     xlim([0,max(Time(end,:))]);
%     ylim([0,20]);
%     LinePlotExplorer();
%     % figure settings for traces
%     set(gca,'Ytick',1:size(reprMaxPos,1)); set(gca,'YtickLabel',identifList(:,1));
%     grid on; xlabel('times, s'); ylabel('trace number');
%     if wavetype==1
%         title('P waves, normalized amplitudes');
%     else
%         title('S waves, normalized amplitudes');
%     end
% end

clear Data;
clear temp identif filenamelist wavetype sampleID hAx;
clear manualPick maxpos;

fprintf('>>> Return from script ''%s'' \n \n', [thisFileName, '.m']);
clear thisFileName thisFolderName;
return

foldpath = 'D:\Geser\Dropbox\IPGG\2016.10.30 - hydrates, analisys\2016.12.04 - wave length analysis\';
save([foldpath,'THF_repr_traces.mat'],'reprTimeP','reprTracesP','reprTimeS','reprTracesS','reprMaxPos','reprDataMIT','reprVpVs','reprQpQs');
