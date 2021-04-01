%% Set paths and directories
clear variables; close all;
% Determine the script name 
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
fprintf('    is running... \n');
addpath(genpath('FuncBase'));
datadir = 'Data\'; % path to data directory
filenamelist = dir(datadir);

%% Parameters for modified Coppens's method
sett = struct( ...
    'windLead', 20,  ... % leading window length (T)
    'windEPS',  30,  ... % EPS operator length (1.5T)
	'coefStab', 0.2, ... % stabilization constant
	'corrType', 2,   ... % 1 - using ER, 2 - using dER, other - no correction
    'skipSmall',  true, ... % true - skip first small amplitude signals
    'checkDistr', true, ... % true - checking distribution of first brakes
    'plotEveryFile', false ... % true - plot all traces, false - plot only error results
    );

% Plotting functions
F2PData = @(x) repmat(1:size(x,2),size(x,1),1) + x./repmat(max(abs(x)),size(x,1),1)/2;
F2PPeaksY = @(x) reshape([-.5+linspace(1,x,x); .5+linspace(1,x,x); nan(1,x)],[],1);
F2PPeaksX = @(x) reshape([reshape(x,1,[]); reshape(x,1,[]); nan(1,length(x))],[],1);

results = struct( ...
    'auto',   NaN(1,size(filenamelist,1)-2), ... % times from autopicking
    'manual', NaN(1,size(filenamelist,1)-2), ... % times from manual picking
    'error',   NaN(1,size(filenamelist,1)-2), ... % absolute values of time differencies
    'meanerr', NaN(1,size(filenamelist,1)-2)  ... % mean of 'error' for corresponding trace
    );

tic
for file = 1:size(filenamelist,1)-2 % first two for '.' and '..'
    fprintf('file: %d from %d \t %s \n',file,size(filenamelist,1)-2,filenamelist(file+2).name);
    
    %% Load data
    load([datadir,filenamelist(file+2).name]); % first two for '.' and '..'
    valuenum = size(Data,1);
    tracenum = size(Data,2);
    
    %% Modified Coppens's method
    % Energy ratio
    DataNorm = Data./repmat(max(Data),size(Data,1),1);
    E1 = filter(ones(1,sett.windLead)',1,DataNorm.^2);
    E2 = cumsum(DataNorm.^2) + sett.coefStab;
    ER = E1./E2;
    clear Data E1 E2;

    % Find first brakes
    AutoPeaks = NaN(tracenum,1);
    BadPeaks  = NaN(tracenum,5);
    EPS = NaN(valuenum,tracenum);
    for trace=1:tracenum
        [AutoPeaks(trace,1), EPS(:,trace)] = fbFindER([(1:size(ER,1))' ER(:,trace)], sett.windEPS, sett.skipSmall, sett.windLead, [(1:size(DataNorm,1))' DataNorm(:,trace)], []);
        if sett.corrType
            AutoPeaks(trace,1) = fbCorrER([(1:size(ER,1))' ER(:,trace)], AutoPeaks(trace,1), sett.windLead, sett.corrType);
        end
    end
    clear trace;

    % Checking distribution of first brakes
    if sett.checkDistr
        for badIter = 1:size(BadPeaks,2)
            deviat = mean(AutoPeaks)-AutoPeaks;
            %deviat = abs(mean(AutoPeaks)-AutoPeaks);
            failpos = find(deviat>min([2.5*std(AutoPeaks),0.7*mean(AutoPeaks)]));
            clear deviat;
            if ~isempty(failpos)
                for i = 1:size(failpos,1)
                    trace = failpos(i,1);
                    BadPeaks(trace,badIter) = AutoPeaks(trace,1);
                    [AutoPeaks(trace,1), EPS(:,trace)] = fbFindER([(1:size(ER,1))' ER(:,trace)], sett.windEPS, sett.skipSmall, sett.windLead, [(1:size(DataNorm,1))' DataNorm(:,trace)], BadPeaks(trace,:));
                    if sett.corrType
                        AutoPeaks(trace,1) = fbCorrER([(1:size(ER,1))' ER(:,trace)], AutoPeaks(trace,1), sett.windLead, sett.corrType);
                    end
                end
                clear trace tracedata failpos deps stdsMat meansMat stds means fbpos i;
            end
        end
        clear badIter failpos;
    end
    
    dEPS = NaN(valuenum,tracenum);
    for i=1:size(EPS,2)
        dEPS(:,i) = circshift(EPS(:,i),-1)-EPS(:,i);
        dEPS(end,i) = NaN;
    end

	%% Save results
	% Error estimation
    err = abs(AutoPeaks-ManualPeaks);
    if tracenum>size(results.auto,1)
        temp = NaN(tracenum-size(results.auto,1),size(results.auto,2));
        results.auto   = [results.auto;   temp];
        results.manual = [results.manual; temp];
        results.error  = [results.error;  temp];
        clear temp;
    end
    results.auto(:,file)    = AutoPeaks;
    results.manual(:,file)  = ManualPeaks;
    results.error(1:size(err,1),file) = err;
    err(isnan(err)) = [];
    errval = mean(err);
    results.meanerr(1,file) = errval;
    
    %% Plotting
    if sett.plotEveryFile
        fig = figure(file);
        set(fig,'units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
        subplot('Position',[0.05 0.07 0.90 0.88]); hold on;
        plot((1:valuenum),F2PData(DataNorm),'k');
        %plot((1:valuenum),F2PData(EPS),'b');
        %plot((1:valuenum),F2PData(dEPS),'g');
        H1 = plot(F2PPeaksX(AutoPeaks),F2PPeaksY(tracenum),'r');
        H2 = plot(F2PPeaksX(ManualPeaks),F2PPeaksY(tracenum),'b');
        legend([H1 H2],'Auto','Manual');
        ylim([0,tracenum+1]); set(gca,'Ytick',1:tracenum);
        axis tight; grid on;
        xlabel('time, discretes'); ylabel('trace, number');
        title([filenamelist(file+2).name,',  mean(abs(\Deltat)) = ',num2str(errval)]);
        clear H1 H2;
        clear err errval;
    end
end
toc

fig = figure;
set(fig,'units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
subplot('Position',[0.05 0.07 0.865 0.88]); hold on;
for file=1:size(results.error,2)
    temp = plot(results.error(:,file),'-o');
    temp.MarkerSize = 9;
end
clear temp;
grid on; leg = legend(num2str((1:size(results.error,2))'));
ylabel('abs(\Deltat)'); xlabel('trace, number');
set(leg,'Position',[0.93 0.33 0.06 0.6]);
clear leg;

fig = figure;
set(fig,'units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
subplot('Position',[0.05 0.07 0.865 0.88]); hold on;
temp = plot(results.meanerr(1,:),'k-o');
temp.MarkerSize = 9;
clear temp;
grid on;
ylabel('mean(abs(\Deltat))'); xlabel('file, number');
xlim([1,size(results.meanerr,2)]); set(gca,'Xtick',1:size(results.meanerr,2));

clear file dne stdNe tracenum valuenum AutoPeaks ManualPeaks;
clear DataNorm ER EPS dER dEPS;
clear F2PData F2PPeaksY F2PPeaksX;

fprintf('>>> Return from script ''%s'' \n \n', [thisFileName, '.m']);
clear thisFolderName thisFileName;
return