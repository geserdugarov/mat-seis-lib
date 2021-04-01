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
datadir = '\\ipgg\project\Проект 14-17\Библиотека экспериментов 2016\Данные\'; % path to directory with data
%datadir = 'data\'; % path to directory with data

%% Settings for V estimation
settV = struct(...
    'firstVal', 0.45, ...    % search first settV.firstVal*max %0.45
    'maxWind',  2.0e-6, ... % s, window after range*max
    'intWind',  0.5e-6, ... % s, window length for interpolatation and finding max
    'polyfitdegr', 3, ...      % degree of interpolation polynom
    'polyfitdisc', 1.0e-9, ... % interpolation discretization
    'noiseFilter', true,                    ... % filter noise or not
    'filtFreqRange', [0.0 0.0 1.5e6 1.6e6], ... % filter before processing ./^^\.
    'useF2Fmean', false ... % use mean correction or curent f2f data for sample
    );

%% Settings for Q estimation
settQ = struct( ...
    'Qestim', true, ... % true - calculate 1/Q
    'dtw', 0.5e-6, ... % s, the length of window for signal smoothing %0.5e-6
    'freqRangeP', [0.15e6, 0.80e6], ... % Hz, frequency range for P wave %[0.2e6, 0.5e6]
    'freqRangeS', [0.15e6, 0.80e6], ... % Hz, frequency range for S wave %[0.2e6, 0.5e6]
... % timewind, s, [shift to start from max, max window from shifted value] 
... % timewind is not required parameter
... 'timewind', [5.0e-6, 15.0e-6], ...
    'numOf0CrossP', 2, ... % number of zero crossing after max position for P %2
    'numOf0CrossS', 2, ... % number of zero crossing after max position for S %4
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
    'saveFigsMAT', false, ... % true - save figure mat-files, false - without saving
    'saveFigsPNG', false, ... % true - save figure png-files, false - without saving
    'saveXLS',  false,  ... % true - save xls file, false - without saving
    'plotXLS',  true,  ... % true - plot in xls file, false - without plotting
    'saveSmpForXLS', true ... % true - save only result png-file near xls-file
    );

%% Mean times for F2F
timesF2F = { ... % [first date after changes, P wave mean max position, S wave mean max position]
    datetime(2016, 4, 7), [4.9482e-06 7.8598e-06], ...
    datetime(2016, 4,13), [4.6887e-06 7.7108e-06], ...
    datetime(2017, 6,10), [5.3453e-06 9.5941e-06]  ...
    };

%% List of sample codes with numbers of bad traces
% 'identificator' [P-wave bad traces] [S-wave bad traces]
% first trace for etalon, second - for aluminum sample
identifList = {
%     'Org' [] []; % plexiglass data
%     'Pes60' [] []; % Pl-Al-La-Cu-St-Cor-Tuf-Pes
    % ! Sb, Sc, Sd and Se have bad discretization, don't do attenuation estimation
%     'Sb' [186 187 422] [186 187 297 422];
%     'Sc' [176] [176];
%     'Sd' [] [];
%     'Se' [13:14 91:92 108:110 112:113] [15 91:93 99 108:110 112:113];
    % ! Sb, Sc, Sd and Se have bad discretization, don't do attenuation estimation
%     'P4128' [] []; % freq. range [0.2 0.8] MHz
%     'P19786' [] []; % freq. range [0.2 0.8] MHz
%     'P4185' [] []; % freq. range [0.2 0.8] MHz
%     'Sh' [15:16 85 317 323] [86 317 323];
%     'Si' [9 408 436] [10 409 436];
%     'Sj' [3:14 138:140] [159]; % from 86 trace heights are equal to 85 measurement
%     'Sn' [4 24:26 99:101 170 240 271:278] [4 103 170 172 241:242 271:278];
%     'So' [11:12 85 155 225 295 365 435 505] [85:86 155:156 225:226 295:296 365:366 436 505:506 519];
%     'Sp' [415 487] [10:11 85 155 225 273 343 487]; % no etalon
%     'Sq' [12 87 157 227 297 367] [12:13 298:299 369];
%     'Sr' [] [];
%     'Ss' [10:11 63 84:85] [12:13 155:156 224:225];
%     'St' [7] [89];
%     'Su' [8 85 242] [11 87:88 157 226:227 242];
%     'Sv' [10 11] [88];
%     'Sw' [4:12 14 15 114:117 182:185 187 238 251:255 258 259 307 308] [6:12 187:189 251 252 257:259];
%     'Sx' [9 83] [9 10 ];
%     'Sy' [10 85 155] [11 64 134 155];
%     'Sac' [8:9 86 134:135 156:157 204 230 275 437] [9:10 136 159 232 297 299 438 439];
%     'Sad' [] [834 1185];
%     'Ice1' [] [1313:1373];
%     'Ice2' [] [388];
%     'SSalta' [] [];
%     'SSaltb' [] [];
%     'SSaltc' [] [];
%     'S002' [23:24 134 343] [24 204 275];
%     'S003' [15:16 97:98] [];
%     'S004' [18 109] [19:20 110];
%     'COAL01' [] [];
%     'claythf' [] [];
%     'SA' [3:8] [3:8];
%     'Sz' [3:8] [3:8];
%     'P4185dry' [3:8] [3:8];
% 	  'claythfs' [3:8] [3:8];
%     'COAL02' [] [];
%     'CH4185' [] [];
 };

%% Saving parameter options in file
if ~exist('outputs','dir')
    mkdir('outputs');
end
save([thisFolderName,'\outputs\settV.mat'],'settV');
save([thisFolderName,'\outputs\settQ.mat'],'settQ');
save([thisFolderName,'\outputs\identifList.mat'],'identifList');
if ~exist('outputs','dir')
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
            searchres = strfind(name,'f2f');
            if isempty(searchres)
                sampname = name;
            else
                etalname = name;
            end
        end
    end
    if ~exist('sampname','var')==1 && exist('etalname','var')==1
        sampname = etalname;
    end
    if exist('sampname','var')==1
        sampdata = load([datadir,sampname]);
    end
    if exist('etalname','var')==1
        etaldata = load([datadir,etalname]);
    end
    % fake f2f for samples without real f2f data
    if ~exist('etalname','var')==1
        fprintf('WARNING: Using fake f2f data for sample: %s\n', sampleID);
        filenamelist2 = dir([thisFolderName,'\data']);
        for i = 1:size(filenamelist2,1)
            name = filenamelist2(i).name;
            searchres = strfind(name,sampleID);
            if ~isempty(searchres)
                etalname = name;
            end
        end
        etaldata  = load([thisFolderName,'\data\',etalname]);
        clear filenamelist2;
    end
    % different load processing for plexiglass
    if strcmp(sampleID,'Org')
        sampdata = load('data\dataOrg.mat'); sampdata = sampdata.dataOrg;
        etaldata = load('data\dataF2F.mat'); etaldata = etaldata.dataF2F;
    end
    if ~exist('sampdata','var')
        error('No sample data was loaded.')
    end
    if ~exist('etaldata','var')
        error('No etalon data was loaded.')
    end

    valuenum = max([size(sampdata.DataP,1),size(sampdata.DataS,1)]);
    tracenum = max([size(sampdata.DataP,2),size(sampdata.DataS,2)])+2;
	clear name searchres i;

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
    % choose f2f trace
    heights = etaldata.DataMIT(9,:);
    pressures = etaldata.DataMIT(5,:);
    [~,posF2F] = max(pressures(heights < 1)); % choose f2f with max pressure
    if isempty(posF2F)
        fprintf('Please check etalon data. Heights for first traces > 1 mm.\n'); 
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
        smplsize = (1.0e-3*sampdata.DataMIT(9,:))'; % sample heights, m
        if settV.useF2Fmean
            % seach right shifts from sample measurement date
            date = datetime(str2double(sampname(1,1:4)),str2double(sampname(1,5:6)),str2double(sampname(1,7:8)));
            dates = [timesF2F{1:2:end}];
            times = [timesF2F{2:2:end}];
            pos = find(date >= dates, 1, 'Last');
            trvltime = maxpos(:,1)-times(wavetype+(pos-1)*2); trvltime(1:2) = []; % traveltimes, s
            clear date dates times pos;
        else
            trvltime = maxpos(:,1)-maxpos(1,1); trvltime(1:2) = []; % traveltimes, s
        end
        Velocities(:,wavetype) = smplsize./trvltime;            % P, S velocities, m/s
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
            end;
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
                if isfield(settQ,'timewind')
                    [~,pos] = min(abs(t0 - ...
                              (maxpos(trace,1)-settQ.timewind(1,1)+settQ.timewind(1,2)) ));
                    temp = cutSignal([Time(:,timeID) Data(:,trace)], ...
                        [maxpos(trace,1)-settQ.timewind(1,1),t0(1,pos)]);
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
                    if any(trace == 2:settQ.plotSkipNum+1:tracenum-1)
                        figure(figID.SigAndSp+(identif-1)*length(fieldnames(figID)));
                        subplot('Position',[(wavetype-1)*0.5+0.05 0.55 0.43 0.41]); hold on;
                        plot(sigCutAft(:,1)-sigCutAft(1,1),sigCutAft(:,2));
                        subplot('Position',[(wavetype-1)*0.5+0.05 0.06 0.43 0.41]); hold on;
                        plot(spAmpAft(:,1),spAmpAft(:,2));
                    end
                end

                switch wavetype
                    case 1
                        freqRange = settQ.freqRangeP;
                    case 2
                        freqRange = settQ.freqRangeS;
                end;
                % calculating Q
                [invQ,invQerr,LnSpOut,linfit] = ...
                    QfromSpRatio(spAmpBef,spAmpAft,freqRange,maxpos(trace,1)-maxpos(1,1),'off','','','');
                if invQ < -1 || invQ > 1
                    invQ = NaN;
                end
                Attenuation(trace-2,wavetype) = invQ;
                % plotting
                if settSavePlot.plotLnSp
                    if any(trace == 3:settQ.plotSkipNum+1:tracenum)
                        figure(figID.lnSp+(identif-1)*length(fieldnames(figID)));
                        subplot('Position',[0.05 (1-wavetype)*0.5+0.56 0.93 0.41]); hold on;
                        plot(freqRange,polyval(linfit,freqRange),'k','LineWidth',1.0);
                        plot(LnSpOut(:,1),LnSpOut(:,2),'-*');
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
    end;
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
    set(FigVQ,'PaperPositionMode','auto');
    saveas(FigVQ,[thisFolderName,'\outputs\VpVsQpQs_',sampleID],'png');
    
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
                       {'process temperature, C'} {'programmed temperature, C'} ...
                       {'temperature inside thermostat, C'} {''} {''} {''} {''}; ...
                       num2cell((1:tracenum-2)') num2cell((sampdata.DataMIT)')];
        else
            if size(sampdata.DataMIT,1)==12
                xlsdata = [{'trace number'} {'time from start, s'} {'temperature inside camera, C'} ...
                           {'pore pressurea at the inlet, mV'} {'pore pressure at the outlet, mV'} ...
                           {'axial pressure, atm'} {'radial pressure, atm'} ...
                           {'-not used-'} {'room temperature, C'} {'sample height, mm'} ...
                           {'process temperature, C'} {'programmed temperature, C'} ...
                           {'temperature inside thermostat, C'}; ...
                           num2cell((1:tracenum-2)') num2cell((sampdata.DataMIT)')];
            else
                xlsdata = [{'trace number'} {'time from start, s'} {'temperature inside camera, C'} ...
                           {'pore pressurea at the inlet, mV'} {'pore pressure at the outlet, mV'} ...
                           {'axial pressure, atm'} {'radial pressure, atm'} ...
                           {'-not used-'} {'room temperature, C'} {'sample height, mm'}; ...
                           num2cell((1:tracenum-2)') num2cell((sampdata.DataMIT)')];
            end;
        end;
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
                if size(sampdata.DataMIT,1)==12
                    eC2.SeriesCollection.NewSeries;
                    eC2.SeriesCollection(2).Value = eS2.Range(['K2:K',num2str(tracenum-1)]);
                    eC2.SeriesCollection(2).XValue = eS.Range(['B2:B',num2str(tracenum-1)]);
                    eC2.SeriesCollection(2).Name = 'T of process';
                    eC2.SeriesCollection(2).Format.Line.ForeColor.RGB = 255*65536 + 0*256 + 0;
                    eC2.SeriesCollection.NewSeries;
                    eC2.SeriesCollection(3).Value = eS2.Range(['L2:L',num2str(tracenum-1)]);
                    eC2.SeriesCollection(3).XValue = eS.Range(['B2:B',num2str(tracenum-1)]);
                    eC2.SeriesCollection(3).Name = 'programmed T';
                    eC2.SeriesCollection(3).Format.Line.ForeColor.RGB = 0*65536 + 0*256 + 255;
                    eC2.SeriesCollection.NewSeries;
                    eC2.SeriesCollection(4).Value = eS2.Range(['M2:M',num2str(tracenum-1)]);
                    eC2.SeriesCollection(4).XValue = eS.Range(['B2:B',num2str(tracenum-1)]);
                    eC2.SeriesCollection(4).Name = 'T inside thermostat';
                    eC2.SeriesCollection(4).Format.Line.ForeColor.RGB = 0*65536 + 0*256 + 0;
                end
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
end;
clear temp identif filenamelist wavetype sampleID hAx;
clear manualPick maxpos;

fprintf('>>> Return from script ''%s'' \n \n', [thisFileName, '.m']);
clear thisFileName thisFolderName;
return