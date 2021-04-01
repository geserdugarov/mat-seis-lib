%% *timePickGather*
% Semi-automatic time picking of microseismic data

%% 
% *Input:* 

% dataset.            - structure containing the following arrays:
%     lenghtUnits     - [string] length units 
%     name            - name of the data set
%     seismicInput    - [noSample, noRec] array of 1C seismic traces 'noSample' samples long  
%                       recorded at 'noRec' receivers
%     time            - [1, noSample] array of time samples
%     timeInput       - [noRec, 1] array of time picks
%     timeRepicked    - [noRec, 1] array of possibly available repicked times
%     timeUnits       - [string] time units ('s' or 'ms')
%     xRec            - [3, noRec] array of receiver locations
%
% processing.         - structure containing the following quantities:
%     contiguous      - [scalar] number of contiguous picks at which such groups of picks are 
%                       rejected
%     maxLag          - [scalar] maximum lag (in timeUnits) for tracking of peaks of troughs
%     rejectAmp       - [scalar] threshold of amplitude of a pick in comparison to 
%                       the maximum amplitude at which the pick is rejected 
%     rejectSNR       - [scalar] threshold of SNR = STA/SLA at which a pick is rejected 
%     timeFlatten     - [scalar] time at which the moveout is flatened
%     windowSTA       - [scalar] half-window around 'timeInput' for computing the average
%                       amplitude of the signal
%     windowLTA       - [scalar] window preceeding 'windowSTA' for computing the average
%                       amplitude of noise
%                       (*) windowSTA/windowLTA is an estimate of signal-to-noise ratio

% fxdecon.            - structure array governing Mauricio Sacchi's f-x deconvolution 
%     freqLow         - [scalar] minimum processed frequency (in Hz) 
%     freqHigh        - [scalar] maximum processed frequency (in Hz) 
%     operatorLength  - [scalar] operator length in samples
%     repetitions     - [scalar] indicating how many times f-x decon is applied to suppress noise 
%                       (default: fxdecon.repetitions = 0)
%     tradeoff        - [scalar] trade-off parameter

% plotting.           - structure containing parameters (mainly from 'plotColorTrace' and 
%                       'plotSeismogram1C')
%     alphaNum        - [scalar] number determining transparancy (see 'plotColorTrace')
%     ampFactor       - [scalar] amplification factor of the traces
%     clipSeismic     - [scalar] from 0 to 1 (default) determining the clip level of seismic
%     figInit         - [scalar] the number of the main figure on which picks are made
%     figPosition     - [2, 4] vector specifying the positions and sizes of figures displaying
%                       seismic  
%     orientation     - [string] indicating the orientation of the time axis:
%                       'vertical' or 'horizontal'
%     peakColor       - [1, 3] color array for plotting the peaks (see 'plotColorTrace')
%     receiverCoord   - [string] indicating 'global' (default) or 'local' coordinate frame 
%                       for plotting the receivers
%     troughColor     - [1, 3] color array for plotting the troughs (see 'plotColorTrace')
%     timeMin         - [scaler] minimum time (in timeUnits) for display 
%     timeMax         - [scaler] maximum time (in timeUnits) for display 
%     time1Color      - [1, 3] color array for plotting the time picks
%     time1MarkerType - [char] marker type for time picks
%     time1LineWidth  - [scalar] line width for time picks
%     time1MarkerSize - [scalar] marker size for time picks
%     time2Color      - [1, 3] color array for plotting the originally repicked times
%     time2MarkerType - [char] marker type for originally repicked times
%     time2LineWidth  - [scalar] line width for originally repicked picks
%     time2MarkerSize - [scalar] marker size for originally repicked times
%     time3MarkerSize - [scalar] marker size for newly repicked times
%     time3MarkerType - [char] marker type for newly repicked times
%     time3LineWidth  - [scalar] line width for newly repicked picks
%     time3MarkerSize - [scalar] marker size for newly repicked times
%     traceLineColor  - [1, 3] color array for lines to plot the traces (see 'plotColorTrace')
%     traceLineWidth  - [scalar] width of lines to plot the traces (see 'plotColorTrace')
%                       (*) if traceLineWidth <= 0, no line is plotted 

%%
% *Output:*

% timeRepicked        - [noRec, 1] array of repicked times 
% indexPT             - [noRec, 1] array of indicators of picked picks (+1) and troughs (-1)

%%
% *Author:* Vladimir Grechka 2013, 2014

%% 
% *Assumptions:*
%
% * the time axis of a seismogram is horizontal
% * the vertical axis is labeled in sequential receiver numbers rather than the receiver coordinates

%%
function [timeRepicked, indexPT] = timePickGather(dataset, processing, fxdecon, plotting)
%% Settings and defaults 
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

markerSizeStepUp = 4;
colorForSeismicImage = 'g';

%% Unpack and check the 'dataset' structure 
noRec = size(dataset.seismicInput, 2);
ampX = NaN(noRec, 1);

time = dataset.time;             
dt = time(2) - time(1);  tmin = time(1);  tmax = time(end);  

sampleInput = round( linMap(dataset.timeInput, tmin, dt, 1, 1) );

if isfield(dataset, 'timeUnits') == 0 || isempty(dataset.timeUnits) == 1
    dataset.timeUnits = 'ms';
end;

if isfield(dataset, 'lengthUnits') == 0 || isempty(dataset.lengthUnits) == 1
    dataset.lengthUnits = 'm';
end;

if isfield(dataset, 'timeInput') == 0 || isempty(dataset.timeInput) == 1
   dataset.timeInput = NaN(noRec, 1);
end;

if isfield(dataset, 'timeRepicked') == 0 || isempty(dataset.timeRepicked) == 1
   dataset.timeRepicked = NaN(noRec, 1);
end;

% Indexes of the output pick/trough indicators
if isfield(dataset, 'indexPT') == 0 || isempty(dataset.indexPT) == 1
    indexPT = NaN(noRec, 1);
else
    indexPT = dataset.indexPT;
end;

% Indexes of existing picks
picks1 = find(isnan(dataset.timeInput) == 0);
picks2 = find(isnan(dataset.timeRepicked) == 0);
picks2T = find(dataset.indexPT == -1);  % troughs

%% Apply defaults to the 'processing' and 'fxdecon' parameters
count = 0;
if isfield(processing, 'timeFlatten') == 0 || isempty(processing.timeFlatten) == 1
    count = count + 1;  processing.timeFlatten = (tmax + tmin)/2 - dt;   
end;

if isfield(processing, 'maxLag') == 0 || isempty(processing.maxLag) == 1
    count = count + 1;  processing.maxLag  = u2u(25, 'ms', dataset.timeUnits);
end;

if isfield(processing, 'windowSTA') == 0 || isempty(processing.windowSTA) == 1
    count = count + 1;  processing.windowSTA = processing.maxLag;
end;

if isfield(processing, 'windowLTA') == 0 || isempty(processing.windowLTA) == 1
    count = count + 1;  processing.windowLTA = 5*processing.windowSTA;
end;

if isfield(processing, 'rejectAmp') == 0 || isempty(processing.rejectAmp) == 1
    count = count + 1;  processing.rejectAmp = 0;
end;

if processing.rejectAmp < 0 || processing.rejectAmp > 1
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Erroneous value of processing.rejectAmp = %g \n', processing.rejectAmp);
    fprintf('>>> It is supposed to be in the interval [0, 1] \n');
    processing.rejectAmp = 0;
    fprintf('>>> Setting processing.rejectAmp = %g -> PAUSE \n', processing.rejectAmp);
    beep on;  beep;  beep off;  pause;
end;

if isfield(processing, 'rejectSNR') == 0 || isempty(processing.rejectSNR) == 1
    count = count + 1;  processing.rejectSNR = 1;
end;

if processing.rejectSNR <= 1
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> The value of processing.rejectSNR = %g \n', processing.rejectSNR);
    fprintf('>>> It is likely too low for the meaningful rejection of picks -> PAUSE \n');
    beep on;  beep;  beep off;  pause;
end;

if isfield(processing, 'contiguous') == 0 || isempty(processing.contiguous) == 1
    count = count + 1;  processing.contiguous = 1;
end;

if isfield(fxdecon, 'repetitions') == 0 || isempty(fxdecon.repetitions) == 1
    count = count + 1;  fxdecon.repetitions = 0;
end;

if count > 0
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    display('>>> Applying defaults to some processing parameters');    
end;

sampleFlatten = round( linMap(processing.timeFlatten, tmin, dt, 1, 1) );

%% Apply defaults for the 'plotting' parameters
lightGray1 = 0.2*[1 1 1];   lightGray2 = 0.85*[1 1 1];    
lightRed  = [1 0.85 0.85];  % lightBlue = [0.5 0.85 1];  
scrsz = get(0, 'ScreenSize');

count = 0;
if isfield(plotting, 'alphaNum') == 0 || isempty(plotting.alphaNum) == 1
    count = count + 1;  plotting.alphaNum = [];
end;

if isfield(plotting, 'ampFactor') == 0 || isempty(plotting.ampFactor) == 1
    count = count + 1;  plotting.ampFactor = 0.9;
end;

if isfield(plotting, 'figInit') == 0 || isempty(plotting.figInit) == 1
    count = count + 1;  plotting.figInit = 1;
end;

if isfield(plotting, 'figPosition') == 0 || isempty(plotting.figPosition) == 1
    count = count + 1;  plotting.figPosition = [50, 50, scrsz(3:4)-100];
end;

if isfield(plotting, 'peakColor') == 0   % isempty(plotting.peakColor) == 1 is legitimate 
    count = count + 1;  plotting.peakColor = lightRed;
end;

if isfield(plotting, 'troughColor') == 0 % isempty(plotting.troughColor) == 1 is legitimate 
    count = count + 1;  plotting.troughColor = lightGray2;  %lightBlue;
end;

if isfield(plotting, 'receiverCoord') == 0 || isempty(plotting.receiverCoord) == 1
    count = count + 1;  plotting.receiverCoord = 'global';  
end;

if isfield(plotting, 'time1Color') == 0 || isempty(plotting.time1Color) == 1
    count = count + 1;  plotting.time1Color = [1 0 0];  % red - 'r'; 
end;

if isfield(plotting, 'time1MarkerType') == 0 || isempty(plotting.time1MarkerType) == 1
    count = count + 1;  plotting.time1MarkerType = 'x'; 
end;

if isfield(plotting, 'time1LineWidth') == 0 || isempty(plotting.time1LineWidth) == 1
    count = count + 1;  plotting.time1LineWidth = 2;
end;

if isfield(plotting, 'time1MarkerSize') == 0 || isempty(plotting.time1MarkerSize) == 1
    count = count + 1;  plotting.time1MarkerSize = 10; 
end;

if isfield(plotting, 'time2Color') == 0 || isempty(plotting.time2Color) == 1
    count = count + 1;  plotting.time2Color = [0 0 1];  % blue - 'b'  
end;

if isfield(plotting, 'time2MarkerType') == 0 || isempty(plotting.time2MarkerType) == 1
    count = count + 1;  plotting.time2MarkerType = '+'; 
end;

if isfield(plotting, 'time2LineWidth') == 0 || isempty(plotting.time2LineWidth) == 1
    count = count + 1;  plotting.time2LineWidth = 2;
end;

if isfield(plotting, 'time2MarkerSize') == 0 || isempty(plotting.time2MarkerSize) == 1
    count = count + 1;  plotting.time2MarkerSize = 12; 
end;

if isfield(plotting, 'time3Color') == 0 || isempty(plotting.time3Color) == 1
    count = count + 1;  plotting.time3Color = [0 1 1];  % cyan - 'c' 
end;

if isfield(plotting, 'time3MarkerType') == 0 || isempty(plotting.time3MarkerType) == 1
    count = count + 1;  plotting.time3MarkerType = 'o'; 
end;

if isfield(plotting, 'time3LineWidth') == 0 || isempty(plotting.time3LineWidth) == 1
    count = count + 1;  plotting.time3LineWidth = 0.5;
end;

if isfield(plotting, 'time3MarkerSize') == 0 || isempty(plotting.time3MarkerSize) == 1
    count = count + 1;  plotting.time3MarkerSize = 3; 
end;

if isfield(plotting, 'traceLineColor') == 0 || isempty(plotting.traceLineColor) == 1
    count = count + 1;  plotting.traceLineColor = lightGray1;  
end;

if isfield(plotting, 'traceLineWidth') == 0 || isempty(plotting.traceLineWidth) == 1
    count = count + 1;  plotting.traceLineWidth = 0.5;  
end;

if isfield(plotting, 'traceOrientation') == 0 || isempty(plotting.traceOrientation) == 1
    count = count + 1;  plotting.traceOrientation = 'horizontal';
elseif strcmp(plotting.traceOrientation, 'vertical') == 1
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    display('>>> The option plotting.traceOrientation = ''vertical'' is not implemented yet');  
    display('>>> Switching to the available option plotting.traceOrientation = ''horizontal'' ');  
    display('>>> PAUSE');  beep on;  beep;  beep off;  pause;  
end;

if count > 0 
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    display('>>> Applying defaults to some plotting parameters');    
end;

%% Plot the map view of receiver layout
figMap = plotting.figInit + 10;
indNotNaN = find(~isnan(dataset.xRec(1,:)));
xRec(1,:) = dataset.xRec(1,:);   
xRec(2,:) = dataset.xRec(2,:);
if strcmp(plotting.receiverCoord, 'global') == 0
    % Plot the map in local coordinates
    xRec(1,:) = dataset.xRec(1,:) - mean(dataset.xRec(1,indNotNaN));  
    xRec(2,:) = dataset.xRec(2,:) - mean(dataset.xRec(2,indNotNaN));
end;
setFigure(figMap, 2, [floor(0.5*scrsz(3:4)), floor(0.43*scrsz(3:4))], ...
    {'Receiver layout', 'Picked peaks (\bullet) and troughs (\o)'}, ...
    ['X (', dataset.lengthUnits, ')'], ...
    ['Y (', dataset.lengthUnits, ')'], [], 14, 'tex'); 
set(gca, 'Ydir', 'normal'); 
if isempty(picks1) == 1
        plot(xRec(1,:), xRec(2,:), '.', 'LineWidth', 0.5, ...       % no original 
        'MarkerEdgeColor', plotting.time1Color, 'MarkerFaceColor', plotting.time1Color);  % picks
else
    plot(xRec(1,picks1), xRec(2,picks1), 'o', ...
        'LineWidth', 0.5, 'MarkerSize', plotting.time3MarkerSize, ...
        'MarkerEdgeColor', plotting.time1Color, 'MarkerFaceColor', plotting.time1Color);  % original
end;
plot(xRec(1,picks2), xRec(2,picks2), 'o', ...
    'LineWidth', 0.5, 'MarkerSize', plotting.time3MarkerSize, ...                     % repicks:
    'MarkerEdgeColor', plotting.time2Color, 'MarkerFaceColor', plotting.time2Color);  % peaks
plot(xRec(1,picks2T), xRec(2,picks2T), 'o', ...
    'LineWidth', 0.5, 'MarkerSize', plotting.time3MarkerSize, ...
    'MarkerEdgeColor', plotting.time2Color, 'MarkerFaceColor', 'w');                  % troughs
axis('equal');  axis('tight');
drawnow;

%% Flatten the moveout at the input time picks and plot the flattened gather
[seismicFlattenAtInputTime, ~] = moveoutCorrect(dataset.seismicInput, ...
    sampleInput, sampleFlatten*ones(size(sampleInput)));

if fxdecon.repetitions > 0
    fprintf('\n>>> Noise suppression with Mauricio Sacchi''s f-x deconvolution. Please wait... \n');    
    for ifx = 1 : floor(fxdecon.repetitions)
        fprintf('>>> Noise-suppression pass %g out of %g \n', ifx, floor(fxdecon.repetitions));
        data0 = zeros(size(seismicFlattenAtInputTime));
        for itrc = 1 : noRec
            trc = seismicFlattenAtInputTime(:,itrc); 
            if max(abs(trc)) > 0
                data0(:,itrc) = trc/max(abs(trc));    % normalize each trace individually
            end;
        end;        

        % Apply Mauricio Sacchi's f-x deconvolution to suppress noise 
        tmp = fx_decon(data0, dt, fxdecon.operatorLength, fxdecon.tradeoff, ...
                                  fxdecon.freqLow,        fxdecon.freqHigh); 
        seismicFlattenAtInputTime = tmp;
    end;
end;

fig = plotting.figInit;   clf(fig);
setFigure(fig, 2, plotting.figPosition(fig,:), ...
    {dataset.name, 'Gather flattened on input (x) time picks; +''s mark original repicked times'}, ...
    ['Time (', dataset.timeUnits, ')'], 'Trace number', [], 14, 'none');  
plotSeismogram1C(fig, time, (1 : noRec)', seismicFlattenAtInputTime, 2, ...
    plotting.ampFactor, plotting.clipSeismic, -1, plotting.timeMin, plotting.timeMax, ...
    processing.timeFlatten*ones(noRec, 1), ...
        plotting.time1Color, plotting.time1MarkerType, ...
        plotting.time1LineWidth, plotting.time1MarkerSize, ...
    dataset.timeRepicked - (dataset.timeInput - processing.timeFlatten), ...
        plotting.time2Color, plotting.time2MarkerType, ...
        plotting.time2LineWidth, plotting.time2MarkerSize, ...
    plotting.traceLineColor, plotting.traceLineWidth, ...  
    plotting.troughColor,    plotting.peakColor, plotting.alphaNum, ...
    plotting.traceOrientation, 'wiggle');   
drawnow;

figImage = plotting.figInit + 1;   clf(figImage);
setFigure(figImage, 2, plotting.figPosition(figImage,:), ...
    {dataset.name, 'Gather flattened on input (x) time picks; green dots mark original repicked times'}, ...
    ['Time (', dataset.timeUnits, ')'], 'Trace number', [], 14, 'none');  
plotSeismogram1C(figImage, time, (1 : noRec)', seismicFlattenAtInputTime, 2, ...
    plotting.ampFactor, plotting.clipSeismic, -1, plotting.timeMin, plotting.timeMax, ...
    [], ... % processing.timeFlatten*ones(noRec, 1), ...
        plotting.time1Color, plotting.time1MarkerType, ...
        plotting.time1LineWidth, plotting.time1MarkerSize, ...
    dataset.timeRepicked - (dataset.timeInput - processing.timeFlatten), ...
        colorForSeismicImage, plotting.time3MarkerType, ...
        plotting.time3LineWidth, plotting.time3MarkerSize, ...
    plotting.traceLineColor, plotting.traceLineWidth, ...  
    plotting.troughColor,    plotting.peakColor, plotting.alphaNum, ...
    plotting.traceOrientation, 'image');  
drawnow;

%% Compute and plot the STA/LTA ratio of amplitudes  
SNR = ratioSTALTA(seismicFlattenAtInputTime, (1 : noRec), sampleFlatten*ones(1, noRec), ...
                  round(processing.windowSTA/dt), round(processing.windowLTA/dt));
dSNR = 0.05*(max(SNR) - min(SNR));

figSNR = plotting.figInit + 4;
setFigure(figSNR, 2, [floor(0.70*scrsz(3)),    plotting.figPosition(1,2), ...
                      floor(0.30*scrsz(3))-10, plotting.figPosition(1,4)], ...
    {'Ratio of amplitudes in short and long', ...
     ['windows: STA = ', num2str(processing.windowSTA), ' ', dataset.timeUnits, ...
             ', LTA = ', num2str(processing.windowLTA), ' ', dataset.timeUnits]}, ...
    'STA-to-LTA amplitude ratio', 'Trace number', [], 14, 'none'); 
plot(SNR, (1 : noRec), 'o-', 'LineWidth', 0.5, 'MarkerSize', 3, ...
    'Color', plotting.time1Color, 'MarkerFaceColor', plotting.time1Color);  
axis([max([0, min(SNR)- dSNR]), max(SNR) + dSNR, 0, noRec + 1]);
    
%% The picking loop

% Settings for the picking loop
pickedTime = NaN(noRec, 1);  pickedTimeSaved   = pickedTime;  
removedPicksSaved = [];      timeRepickedSaved = [];
pickedTraceNumber = [1, noRec];

clear tmp
noMax = round(processing.maxLag/dt);  % max allowable time difference (in samples)

while pickedTraceNumber(1) ~= 0 
    
    %% Making picks
    errorFlag = 1;  
    while errorFlag == 1 || isempty(pickedTraceNumber) == 1  
        try 
            errorFlag = 0;
            fprintf('\n>>> Please type a range of trace numbers [traceBeg, traceEnd] \n');
            pickedTraceNumber = round( input( ...
                '    to manually pick trace ''traceBeg'' or enter 0 (zero) to exit -> ') );         
        catch err; 
            errorFlag = 1;  
            fprintf('>>> %s \n', err.message);
            display('>>> Please correct the input'); 
        end; 
        
        if isempty(pickedTraceNumber) == 1  
            errorFlag = 1;  
            display('>>> No trace number was entered'); 
            display(['>>> Possible repeated messages, caused by accumulation of ''ENTER''-type ', ...
                     'actions during event picking and removal, are to be disregarded']); 
        elseif min(pickedTraceNumber) < 0 || max(pickedTraceNumber) > noRec
            errorFlag = 1;  
            fprintf('>>> Erroneous trace number [%g, %g] -> it should be between 1 and noRec = %g \n', ...
                    pickedTraceNumber, noRec);
            display('>>> Please correct the input'); 
        end;
        
        if length(pickedTraceNumber) > 2
            errorFlag = 1;  
            display('>>> The length of input should be 1 or 2'); 
        end;
        
        if length(pickedTraceNumber) == 2  &&  pickedTraceNumber(1) > pickedTraceNumber(2)
            errorFlag = 1;  
            display('>>> The first input trace number should be smaller than or equal to the second'); 
        end;
        
        if length(pickedTraceNumber) == 1  
            pickedTraceNumber(2) = pickedTraceNumber(1); 
            errorFlag = 0;  
            if pickedTraceNumber(1) ~= 0 
                display('>>> WARNING: Single-number input - ''traceEnd'' is set equal to ''traceBeg'' '); 
            end;
        end;
    end;
    
    if pickedTraceNumber(1) == 0;  break;  end;
    
    fprintf('>>> Please pick time at trace %g \n', pickedTraceNumber(1));
    figure(fig);
    [timX, traceX] = ginput(1); 
    pickedTime(pickedTraceNumber(1)) = timX;            

    if pickedTraceNumber(1) - 1 <= round(traceX) && round(traceX) <= pickedTraceNumber(1) + 1 
        ampX(pickedTraceNumber(1), 1) = ...
            seismicFlattenAtInputTime(round(linMap(timX, tmin, dt, 1, 1)), round(traceX)); 
        if ampX(pickedTraceNumber(1), 1) < 0
            indexPT(pickedTraceNumber(1), 1) = 1;    % peak
        else
            indexPT(pickedTraceNumber(1), 1) = -1;   % trough
        end;

        % Automatically track the selected phase through the gather from 
        % 'pickedTraceNumber(1)' to 'pickedTraceNumber(2)'
        noPickAccept = 0;  noPickReject = 0;
        pickedSampleLast = round(linMap(pickedTime(pickedTraceNumber(1)), tmin, dt, 1, 1));
        for jrec = pickedTraceNumber(1) + 1 : pickedTraceNumber(2)
            pickedSampleBeg  = pickedSampleLast - noMax;
            pickedSampleEnd  = pickedSampleLast + noMax;
            tmp = seismicFlattenAtInputTime(pickedSampleBeg : pickedSampleEnd, jrec);

            % This block maintains polarity of the manually made pick 
            if indexPT(pickedTraceNumber(1), 1) == 1        % peak
                sampleMaxMin = find(tmp == min(tmp));
            elseif indexPT(pickedTraceNumber(1), 1) == -1   % trough
                sampleMaxMin = find(tmp == max(tmp));
            end;
            pickedSample = pickedSampleLast - (noMax + 1) + sampleMaxMin(1);

            % Compute the STA/SLA ratio
            sampleCenter = NaN(noRec,1);  sampleCenter(jrec,1) = pickedSample;
            currentSNR = ratioSTALTA(seismicFlattenAtInputTime, jrec, sampleCenter, ...
                                     round(processing.windowSTA/dt), ...
                                     round(processing.windowLTA/dt));      

            if abs(tmp(sampleMaxMin(1)))/abs(ampX(pickedTraceNumber(1), 1)) > processing.rejectAmp && ...
               currentSNR(jrec) > processing.rejectSNR
                % Accepted pick
                pickedSampleLast = pickedSample;                            
                pickedTime(jrec) = linMap(pickedSample, 1, 1, tmin, dt);
                indexPT(jrec, 1) = indexPT(pickedTraceNumber(1), 1);
                fprintf(['>>> Receiver = %4.0f: SNR = %5.3f, amplitude ratio = %5.3f -> ', ...
                         'The pick is accepted, time = %g (%s) \n'], jrec, ...
                    currentSNR(jrec), abs(tmp(sampleMaxMin(1)))/abs(ampX(pickedTraceNumber(1), 1)), ...
                    linMap(pickedSample, 1, 1, tmin, dt), dataset.timeUnits);
                noPickAccept = noPickAccept + 1;
            else
                % Rejected pick
                pickedTime(jrec) = NaN;
                indexPT(jrec, 1) = NaN;
                fprintf(['>>> Receiver = %4.0f: SNR = %5.3f, amplitude ratio = %5.3f -> ', ...
                         'The pick is rejected \n'], jrec, ...
                    currentSNR(jrec), abs(tmp(sampleMaxMin(1)))/abs(ampX(pickedTraceNumber(1), 1)));
                noPickReject = noPickReject + 1;
            end;        
        end;

        % Remove groups containing 'processing.contiguous' or fewer contiguous picks
        jcontig = 0;
        for jrec = pickedTraceNumber(1) + 1 : pickedTraceNumber(2)
            if 0 < jcontig && jcontig <= processing.contiguous && isnan(pickedTime(jrec)) == 1
                pickedTime(jrec-jcontig : jrec-1) = NaN;
                indexPT(jrec-jcontig : jrec-1, 1) = NaN;
            end;
            if isnan(pickedTime(jrec)) == 1
                jcontig = 0;
            else
                jcontig = jcontig + 1;  % accumulate the number of contiguous picks
            end;
        end;
        fprintf('\n>>> Statistics: %g accepted picks, %g rejected picks \n \n', ...
                noPickAccept, noPickReject);

        % Plot the newly made picks
        figure(fig);
        plot(pickedTimeSaved, 1:noRec, plotting.time3MarkerType, ...
            'LineWidth', plotting.time3LineWidth, ...
            'MarkerSize', plotting.time3MarkerSize + markerSizeStepUp, ...
            'MarkerEdgeColor', 'w', 'MarkerFaceColor', 'w');
        plot(pickedTime,      1:noRec, plotting.time3MarkerType, ...
            'LineWidth', plotting.time3LineWidth, ...
            'MarkerSize', plotting.time3MarkerSize + markerSizeStepUp, ...
            'MarkerEdgeColor', 'w', 'MarkerFaceColor', plotting.time3Color);

        % Update the image
        figure(figImage);
        plot(pickedTimeSaved, 1:noRec, plotting.time3MarkerType, ...
            'LineWidth', plotting.time3LineWidth, ...
            'MarkerSize', plotting.time3MarkerSize, ...
            'MarkerEdgeColor', 'w', 'MarkerFaceColor', 'w');
        plot(pickedTime,      1:noRec, plotting.time3MarkerType, ...
            'LineWidth', plotting.time3LineWidth, ...
            'MarkerSize', plotting.time3MarkerSize, ...
            'MarkerEdgeColor', plotting.time3Color, 'MarkerFaceColor', plotting.time3Color);

        % Update the map
        figure(figMap);
        inopick = find( isnan(indexPT(pickedTraceNumber(1) : pickedTraceNumber(2))) == 1 );
        ipeak   = find( indexPT(pickedTraceNumber(1) : pickedTraceNumber(2)) == +1 );
        itrough = find( indexPT(pickedTraceNumber(1) : pickedTraceNumber(2)) == -1 );
        plot(xRec(1, pickedTraceNumber(1) + inopick - 1), ...  % original picks
             xRec(2, pickedTraceNumber(1) + inopick - 1), ... 
             'o', 'LineWidth', plotting.time3LineWidth, ...
             'MarkerSize', plotting.time3MarkerSize, ...
             'MarkerEdgeColor', plotting.time1Color, 'MarkerFaceColor', plotting.time1Color);  
        plot(xRec(1, pickedTraceNumber(1) + ipeak - 1), ...    % newly made peaks
             xRec(2, pickedTraceNumber(1) + ipeak - 1), ...
             'o', 'LineWidth', plotting.time3LineWidth, ...
             'MarkerSize', plotting.time3MarkerSize, ...
             'MarkerEdgeColor', plotting.time3Color, ...
             'MarkerFaceColor', plotting.time3Color);
        plot(xRec(1, pickedTraceNumber(1) + itrough - 1), ...  % newly made troughs
             xRec(2, pickedTraceNumber(1) + itrough - 1), ...
             'o', 'LineWidth', plotting.time3LineWidth, ...
             'MarkerSize', plotting.time3MarkerSize, ...
             'MarkerEdgeColor', plotting.time3Color, 'MarkerFaceColor', 'w');

        %% Interactive pick removal
        fprintf('>>> Interactive removal of mispicks \n');
        mouseButton = 1;
        while mouseButton == 1
            fprintf(['>>> When the cross-hair is displayed, the pick-removal ', ...
                     'process is ON. In this mode, \n']);
            fprintf(['    (*) clicking on a trace with the LEFT mouse button ', ...
                     'removes pick from it \n']);
            fprintf(['    (*) ''ENTER'' accepts the removed picks, if any, ', ...
                     'and exits from the current stage of the pick removal \n']);
            fprintf(['    (*) the RIGHT mouse-button click followed by ''ENTER'' ', ...
                     'terminate the pick-removal process \n']);

            removedPicks = [];  
            figure(fig);
            [~, traceRemoved, tmpMouseButton] = ginput;  % manual removal of picks
            % Note: Picks can be removed in any order. Repeated clicking on a trace or 
            %       clicking on a trace that has no pick is fine but clicking outside 
            %       the trace range causes error.

            % Analyze the clicks
            findRemoved = find(tmpMouseButton == 1);
            findStop = find(tmpMouseButton == 3, 1);

            if isempty(findRemoved) == 0
                fprintf('>>> Removing picks from trace(s) \n');
                disp(round(traceRemoved(findRemoved)'));

                removedPicks(1 : length(findRemoved)) = round(traceRemoved(findRemoved));

                figure(fig);
                plot(pickedTime(removedPicks), removedPicks, plotting.time3MarkerType, ...
                    'LineWidth', plotting.time3LineWidth, ...
                    'MarkerSize', plotting.time3MarkerSize + markerSizeStepUp, ...
                    'MarkerEdgeColor', 'w', 'MarkerFaceColor', 'w');
                figure(figImage);
                plot(pickedTime(removedPicks), removedPicks, plotting.time3MarkerType, ...
                    'LineWidth', plotting.time3LineWidth, ...
                    'MarkerSize', plotting.time3MarkerSize, ...
                    'MarkerEdgeColor', 'w', 'MarkerFaceColor', 'w');
                pickedTime(removedPicks) = NaN;                     

                figure(fig);
                plot(dataset.timeRepicked(removedPicks) - ...
                    (dataset.timeInput(removedPicks) - processing.timeFlatten), ...
                    removedPicks, ...
                    plotting.time2MarkerType, 'LineWidth', plotting.time2LineWidth, ...
                    'MarkerSize', plotting.time2MarkerSize + 1, ...
                    'MarkerEdgeColor', 'w', 'MarkerFaceColor', 'w');
                figure(figImage);
                plot(dataset.timeRepicked(removedPicks) - ...
                    (dataset.timeInput(removedPicks) - processing.timeFlatten), ...
                    removedPicks, ...
                    plotting.time3MarkerType, 'LineWidth', plotting.time3LineWidth, ...
                    'MarkerSize', plotting.time3MarkerSize, ...
                    'MarkerEdgeColor', 'w', 'MarkerFaceColor', 'w');

                removedPicksSaved = cat(2, removedPicksSaved, removedPicks);
                timeRepickedSaved = cat(2, timeRepickedSaved, ...
                                           dataset.timeRepicked(removedPicks)');
                dataset.timeRepicked(removedPicks) = NaN;
                indexPT(removedPicks, 1) = NaN;

                % Redraw the initial picks
                figure(fig);
                plot(processing.timeFlatten*ones(size(removedPicks)), removedPicks, ...
                    plotting.time1MarkerType, 'LineWidth', plotting.time1LineWidth, ...
                    'MarkerSize', plotting.time1MarkerSize, ...
                    'MarkerEdgeColor', plotting.time1Color, ...
                    'MarkerFaceColor', plotting.time1Color);

                % Update the map
                figure(figMap);
                plot(xRec(1,removedPicks), xRec(2,removedPicks), ...
                    'o', 'LineWidth', plotting.time3LineWidth, ...
                    'MarkerSize', plotting.time3MarkerSize, ...
                    'MarkerEdgeColor', 'w', 'MarkerFaceColor', 'w');
                plot(xRec(1,removedPicks), xRec(2,removedPicks), ...
                    'o', 'LineWidth', 0.5, 'MarkerSize', 3, ...
                    'MarkerEdgeColor', plotting.time1Color, ...
                    'MarkerFaceColor', plotting.time1Color);
            end;    

            if isempty(findStop) == 1
                mouseButton = 1;
                fprintf('>>> Control is returned to Figure %g \n', fig);
                fprintf(['>>> To return to the pick-removal mode, please click on ', ...
                         'the Command Window and then ''ENTER'' \n']);  pause;
            else
                 mouseButton = 3;
                 fprintf('>>> Exiting pick-removal process \n');
            end;
        end;  % of the pick-removal while-loop

        %% Interactive interpolation of picks 
        fprintf('\n>>> Interactive interpolation of picks \n \n');

        traceBeg = 1;
        while traceBeg ~= 0
            errorFlag = 1;  % try to catch input errors
            while errorFlag == 1 
                errorFlag = 0;
                traceBeg = input(['>>> Please input 1 (one) to begin linear ', ...
                                  'interpolation of picks or 0 (zero) to exit -> ']); 
                if isempty(traceBeg) == 1 || isnumeric(traceBeg(1)) == 0 || ...
                   (traceBeg(1) ~= 0 && traceBeg(1) ~= 1)
                    errorFlag = 1;  
                    fprintf(['>>> The first element of input = %g, ', ...
                             'whereas its acceptable values are 0 or 1 \n'], traceBeg(1)); 
                    fprintf('>>> Please correct the input \n \n'); 
                end; 
            end;  % of catching an error 
            
            if traceBeg == 0;  break;  end;

            fprintf(['>>> When the cross-hair is displayed, left-click sets the nodes ', ...
                     'for linear interpolation and ''ENTER'' terminates the process \n']);
            fprintf('>>> The interpolation starts at the first click and ends at the last \n');
            fprintf(['    (*) Please wait until the the cross-hair is properly displayed ', ...
                     'because clicks might not be accepted otherwise \n']);
            fprintf(['    (*) Clicking on another window when in interpolation mode often ', ...
                     'creates erroneous picks \n \n']);
            
            clear tmpTrace tmpPick tmpTime tmpPT
            figure(figImage);
            [tmpPick0, tmpTrace0] = ginput;             
            [tmpTrace, tmpTraceOrder] = sort(tmpTrace0, 'ascend');  % sort the picks in 
            tmpPick = tmpPick0(tmpTraceOrder);                      % ascending order
            traceBeg = round(tmpTrace(1));
            traceEnd = round(tmpTrace(end));
            tmpTime = interp1(tmpTrace, tmpPick, (traceBeg : traceEnd), 'linear', 'extrap');
            fprintf('>>> Linear interpolation from picks \n');
            disp(tmpPick');
            fprintf('    at traces \n');
            disp(round(tmpTrace'));
                
%             traceEnd = round( input(['>>> Please input number of the last  trace ', ...
%                                      'for interpolation -> ']) ); 
%             clear tmpRecs tmpTime tmpPT
%             tmpRecs = (traceBeg : traceEnd)';
%             tmpTime = pickedTime(tmpRecs);
%             tmpRecs(isnan(tmpTime)) = [];  % remove NaNs
%             tmpTime(isnan(tmpTime)) = [];
%             tmpPol = polyfit(tmpRecs, tmpTime, 2);  % quadratic smoother/infiller             
%             clear tmpTime
%             tmpTime = polyval(tmpPol, (traceBeg : traceEnd)');

            for itrace = traceBeg : traceEnd 
%                itrace
%                round(linMap(tmpTime(itrace - traceBeg + 1), tmin, dt, 1, 1))
                ampX(itrace - traceBeg + 1, 1) = seismicFlattenAtInputTime( ...
                    round(linMap(tmpTime(itrace - traceBeg + 1), tmin, dt, 1, 1)), round(itrace)); 
                if ampX(itrace - traceBeg + 1, 1) < 0
                    indexPT(itrace, 1) = 1;    % peak
                else
                    indexPT(itrace, 1) = -1;   % trough
                end;
%                indexPT(itrace, 1)
            end;
        
            % Updata the wiggle trace display
            figure(fig);
            plot(pickedTime(traceBeg : traceEnd), (traceBeg : traceEnd), plotting.time3MarkerType, ...
                'LineWidth', plotting.time3LineWidth, ...
                'MarkerSize', plotting.time3MarkerSize + markerSizeStepUp, ...
                'MarkerEdgeColor', 'w', 'MarkerFaceColor', 'w');
            plot(tmpTime, (traceBeg : traceEnd), plotting.time3MarkerType, ...
                'LineWidth', plotting.time3LineWidth, ...
                'MarkerSize', plotting.time3MarkerSize + markerSizeStepUp, ...
                'MarkerEdgeColor', 'w', 'MarkerFaceColor', plotting.time3Color);

            % Update the image
            figure(figImage);
            plot(pickedTime(traceBeg : traceEnd), (traceBeg : traceEnd), plotting.time3MarkerType, ...
                'LineWidth', plotting.time3LineWidth, ...
                'MarkerSize', plotting.time3MarkerSize, ...
                'MarkerEdgeColor', 'w', 'MarkerFaceColor', 'w');
            plot(tmpTime, (traceBeg : traceEnd), plotting.time3MarkerType, ...
                'LineWidth', plotting.time3LineWidth, ...
                'MarkerSize', plotting.time3MarkerSize, ...
                'MarkerEdgeColor', plotting.time3Color, 'MarkerFaceColor', plotting.time3Color);

            % Update the map
            figure(figMap);
            for itrace = traceBeg : traceEnd 
                if indexPT(itrace, 1) == 1
                    plot(xRec(1, traceBeg : traceEnd), ...    % interpolated  peaks
                         xRec(2, traceBeg : traceEnd), ...
                         'o', 'LineWidth', plotting.time3LineWidth, ...
                         'MarkerSize', plotting.time3MarkerSize, ...
                         'MarkerEdgeColor', plotting.time3Color, ...
                         'MarkerFaceColor', plotting.time3Color);
                else
                    plot(xRec(1, traceBeg : traceEnd), ...    % interpolated troughs
                         xRec(2, traceBeg : traceEnd), ...
                         'o', 'LineWidth', plotting.time3LineWidth, ...
                         'MarkerSize', plotting.time3MarkerSize, ...
                         'MarkerEdgeColor', plotting.time3Color, 'MarkerFaceColor', 'w');
                end;
            end;
            pickedTime(traceBeg : traceEnd) = tmpTime;
        end;  % of interpolation while-loop
        fprintf('>>> Exiting interactive interpolation of picks \n');
        pickedTimeSaved = pickedTime;       

    else
        % The pick is made on incorrect trace
        fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
        fprintf(['>>> The pick is made on trace %g rather than on the selected ', ...
                'trace %g \n'], ...
                round(traceX), pickedTraceNumber(1));
        fprintf('>>> Please repick \n');  
        beep on;  beep;  beep off;
    end;
    
    if pickedTraceNumber(1) ~= 0 
        fprintf('\n>>> Last manual pick was made on trace number %g \n', pickedTraceNumber(1));
    else
        fprintf('\n>>> Exiting time picking \n'); 
    end;
    
end;  % of the main while-loop

%% Compute samples corresponding to the original and new repicked times
pickedSample = NaN(noRec, 1);
for irec = 1:noRec
    if isnan(pickedTime(irec)) == 0
        pickedSample(irec) = round( linMap(pickedTime(irec), tmin, dt, 1, 1) );
        dataset.timeRepicked(irec) = NaN;
    else
        if isnan(dataset.timeRepicked(irec)) == 0
            pickedSample(irec) = round( linMap( ...
            dataset.timeRepicked(irec) - (dataset.timeInput(irec) - processing.timeFlatten), ...
            tmin, dt, 1, 1) );
        end;
    end;
end;

%% Moveout-correct and plot the repicked gather
[seismicFlattenAtPickedTime, ~] = moveoutCorrect(seismicFlattenAtInputTime, ...
    pickedSample, sampleFlatten*ones(size(pickedSample)));    

fig = fig + 20;
setFigure(fig, 2, plotting.figPosition(1,:), ...
    {dataset.name, 'Gather flattened on repicked times (o and +)'}, ...
    ['Time (', dataset.timeUnits, ')'], 'Trace number', [], 14, 'none');  
plotSeismogram1C(fig, time, (1 : noRec)', seismicFlattenAtPickedTime, 2, ...
    plotting.ampFactor, plotting.clipSeismic, -1, plotting.timeMin, plotting.timeMax, ...
    processing.timeFlatten*ones(noRec, 1) + 0*dataset.timeRepicked, ...
        plotting.time2Color, plotting.time2MarkerType, ...
        plotting.time2LineWidth, plotting.time2MarkerSize, ...
    2*processing.timeFlatten - pickedTime, ...
        plotting.time1Color, plotting.time1MarkerType, ...
        plotting.time1LineWidth, plotting.time1MarkerSize, ...
    plotting.traceLineColor, plotting.traceLineWidth, ...  
    plotting.troughColor,    plotting.peakColor, plotting.alphaNum, ...
    plotting.traceOrientation, 'wiggle');   
plot(processing.timeFlatten*ones(noRec, 1) + 0*pickedTime, 1:noRec, plotting.time3MarkerType, ...
    'LineWidth',  plotting.time3LineWidth,  'MarkerSize', plotting.time3MarkerSize, ...
    'MarkerEdgeColor', plotting.time3Color, 'MarkerFaceColor', plotting.time3Color);

fig = fig + 21;
setFigure(fig, 2, plotting.figPosition(2,:), ...
    {dataset.name, 'Gather flattened on repicked times (o and +)'}, ...
    ['Time (', dataset.timeUnits, ')'], 'Trace number', [], 14, 'none');  
plotSeismogram1C(fig, time, (1 : noRec)', seismicFlattenAtPickedTime, 2, ...
    plotting.ampFactor, plotting.clipSeismic, -1, plotting.timeMin, plotting.timeMax, ...
    [], ...
        plotting.time2Color, plotting.time2MarkerType, ...
        plotting.time2LineWidth, plotting.time2MarkerSize, ...
    [], ...
        plotting.time1Color, plotting.time1MarkerType, ...
        plotting.time1LineWidth, plotting.time1MarkerSize, ...
    plotting.traceLineColor, plotting.traceLineWidth, ...  
    plotting.troughColor,    plotting.peakColor, plotting.alphaNum, ...
    plotting.traceOrientation, 'image');   
%plot(processing.timeFlatten*ones(noRec, 1) + 0*pickedTime, 1:noRec, plotting.time3MarkerType, ...
%    'LineWidth',  plotting.time3LineWidth,  'MarkerSize', plotting.time3MarkerSize, ...
%    'MarkerEdgeColor', colorForSeismicImage, 'MarkerFaceColor', colorForSeismicImage);

%% Compute the new picks 
timeRepicked = pickedTime + (dataset.timeInput - processing.timeFlatten);

% Fill array 'timeRepicked' with the originally repicked times
for irec = 1:noRec
    if isnan(timeRepicked(irec)) == 1
        timeRepicked(irec) = dataset.timeRepicked(irec);
    end;
end;

% %% Display the new picks on original data 
% figure(plotting.figInit);
% plot(timeRepicked, 1:noRec, plotting.time3MarkerType, 'LineWidth', plotting.time3LineWidth, ...
%     'MarkerSize', plotting.time3MarkerSize, ...
%     'MarkerEdgeColor', 'w', 'MarkerFaceColor', plotting.time3Color);
% 
% plot(timeRepickedSaved, removedPicksSaved, ...
%     plotting.time2MarkerType, 'LineWidth', plotting.time2LineWidth, ...
%     'MarkerSize', plotting.time2MarkerSize + 1, 'MarkerEdgeColor', 'w', 'MarkerFaceColor', 'w');

end  % of the function
