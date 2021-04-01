%% *plotSeismogram3C*
% Plot three-component seismogram

%%
% *Input:*

% fig                - [scalar] number of a new figure
% time               - [1, noTime] array of time samples
%                      (*) if array time is empty, the traces are ploted in samples
% recCoor            - [noRec, 1] - array of receiver coordinates 
%                      (*) if array recCoor is empty, the seismogram is labeled in trace numbers
% traceDataT         - [3*noRec, noTime] matrix containing traces
% flagNorm           - [scalar] normalization flag
%                      flagNorm = 0 - no normalization (default)
%                      flagNorm = 1 - normalize to the maximim of the entire seismogram
%                      flagNorm = 2 - normalize each trace component individually
%                      flagNorm = 3 - normalize each three-component trace 
% ampFactor          - [scalar] amplification factor of the traces
%                      (*) if variable ampFactor is empty, ampFactor is defaulted to 1
% polarity           - [scalar] equal to +1 or -1 to indicate the peak/trough switch
% tmin               - [scalar] minimum time for plotting
%                      (*) if variable tmin is empty, tmin = min(time)
% tmax               - [scalar] maximum time for plotting
%                      (*) if variable tmax is empty, tmax = max(time)
% timePick           - [noRec, :] array of picked times
% timePickColor      - [1, 3] color array for plotting the time picks
% timePickMarkerType - [char] marker type for time picks
% timePickLineWidth  - [scalar] line width for time picks
% timePickMarkerSize - [scalar] marker size for time picks
% timeComp           - [noRec, :] array of computed times
% timeCompColor      - [1, 3] color array for plotting the computed times
% timeCompMarkerType - [char] marker type for computed times
% timeCompLineWidth  - [scalar] line width for computed picks
% timeCompMarkerSize - [scalar] marker size for computed times
% traceLineWidth     - [scalar] width of lines to plot the traces (see 'plotColorTrace')
%                      (*) if traceLineWidth <= 0, no line is plotted 

%%
% *Output:* none

%%
% *Author:* Vladimir Grechka 2012 - 2014

%%
function plotSeismogram3C(fig, time, recCoor, traceDataT, flagNorm, ampFactor, polarity, tmin, tmax, ...
    timePick, timePickColor, timePickMarkerType, timePickLineWidth, timePickMarkerSize, ...
    timeComp, timeCompColor, timeCompMarkerType, timeCompLineWidth, timeCompMarkerSize, traceLineWidth)
%% Settings 
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

traceData = traceDataT';
noRec = floor(size(traceData, 1)/3);
noTime = size(traceData, 2);
clipLevel = 0.95;

%% Checks and defaults
if isempty(time) == 1
    time(1,:) = (1:noTime);                         % default for variable 'time'
end;
if size(time, 2) ~= noTime
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> size(time) = [%g, %g],  size(traceData) = [%g, %g] \n', ...
            size(time), size(traceData));    
    fprintf('>>> Inconsistent sizes of arrays ''time'' and ''traceData'' \n \n');    
      error('>>> STOP');
end;

if isempty(recCoor) == 1
    recCoor(:,1) = (1:noRec)';                      % default for variable 'recCoor'
end;
if size(recCoor, 1) ~= noRec
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> size(recCoor) = [%g, %g],  size(traceData) = [%g, %g] \n', ...
            size(recCoor), size(traceData));    
    fprintf('>>> Inconsistent sizes of arrays ''recCoor'' and ''traceData'' \n \n');    
      error('>>> STOP');
end;

if isequal(size(traceData, 1)/3, noRec) == 0
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> The number of traces in 3C seismogram (%g) is not divisible by 3 \n \n', ...
            size(traceData, 1));
      error('>>> STOP');
end; 

if isempty(flagNorm) == 1  ||  sum(flagNorm == [0, 1, 2, 3]) == 0
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf(['>>> WARNING: Variable ''flagNorm'' is defaulted to 0 because it is either ', ...
             'empty or its value [= %g] is unsupported \n'], flagNorm);    
    flagNorm = 0;
end;

if isempty(tmin) == 1   
    tmin = time(1,1);                               % default for variable 'tmin'
end;       

if isempty(tmax) == 1   
    tmax = time(1,noTime);                          % default for variable 'tmax'
end;  

if isempty(flagNorm) == 1   
    flagNorm = 0;                                   % default for variable 'flagNorm'
end;  

if isempty(ampFactor) == 1   
    ampFactor = 1;                                  % default for variable 'ampFactor'
end;  

if isempty(timePickColor) == 1
    timePickColor = 'b';                            % default for variable 'timePickColor'
end;

if isempty(timePickMarkerType) == 1
    timePickMarkerType = '+';                       % default for variable 'timePickMarkerType'
end;

if isempty(timePickLineWidth) == 1
    timePickLineWidth = 0.5;                        % default for variable 'timePickLineWidth'
end;

if isempty(timePickMarkerSize) == 1
    timePickMarkerSize = 10;                        % default for variable 'timePickMarkerSize'
end;

if isempty(timeCompColor) == 1
    timeCompColor = 'r';                            % default for variable 'timeCompColor'
end;

if isempty(timeCompMarkerType) == 1
    timeCompMarkerType = 'x';                       % default for variable 'timeCompMarkerType'
end;

if isempty(timeCompLineWidth) == 1
    timeCompLineWidth = 0.5;                        % default for variable 'timeCompLineWidth'
end;

if isempty(timeCompMarkerSize) == 1
    timeCompMarkerSize = 10;                        % default for variable 'timeCompMarkerSize'
end;

%% Set up plotting
figure(fig);
grid on;   set(gca, 'YGrid', 'off');   hold on; 

xlim([tmin, tmax]);
if numel(recCoor) > 1 
    dCoor = abs(recCoor(2) - recCoor(1));
else
    dCoor = 1;
end;
ylim([min(recCoor) - 0.999*dCoor, max(recCoor) + 0.999*dCoor]);

if polarity == 1
    set(gca, 'YDir', 'normal');    
elseif polarity == -1;
    set(gca, 'YDir', 'reverse');      
else
    set(gca, 'YDir', 'normal');    
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> WARNING: Incorrect value of polarity = %g \n', polarity);    
    fprintf('>>> Variable ''YDir'' is defaulted to to ''normal'' \n \n');    
end;

%% Normalize traces
if flagNorm == 1
    maxTrace = max(max(abs(traceData)));
    if maxTrace > 0
        traceData = traceData/maxTrace;             % normalize the entire seismogram
    end;
end;

%% Plot the seismogram
red = [1, 0, 0];    gray = 0.5*[1, 1, 1];    blue = [0, 0, 1];
traceLineColor = [red; gray; blue];
updw = 0.1*[-1, 0, 1];
for itrc = 1:noRec
    trc3C = traceData(3*(itrc - 1) + 1 : 3*itrc, :);
    if flagNorm == 3
        max3Ctrace = max(max(abs(trc3C)));
        if max3Ctrace > 0
            trc3C = trc3C/max3Ctrace;               % normalize each three-component trace 
        end;
    end;
        
    for icomp = 1:3
        trc = trc3C(icomp, :); 
        if flagNorm == 2  &&  max(abs(trc)) > 0
            trc = trc/max(abs(trc));                % normalize each trace component individually
        end;
        
        trc = ampFactor*trc;                        % apply the amplitude factor
        trc(trc >  clipLevel) =  clipLevel;         % clip the traces
        trc(trc < -clipLevel) = -clipLevel;
    
        plotColorTrace(time, recCoor(itrc) + dCoor*(updw(icomp) + trc), ...
                       traceLineColor(icomp,:), traceLineWidth, [], [], []);
    end;

    % Plot the computed and picked times
    if isempty(timeComp) == 0
        for iwave = 1:size(timeComp, 2)
            if isempty(timeComp(itrc, iwave)) == 0
                plot(timeComp(itrc, iwave), recCoor(itrc), ...
                     timeCompMarkerType, 'LineWidth', timeCompLineWidth, ...
                     'MarkerSize', timeCompMarkerSize, 'Color', timeCompColor);   
            end;       
        end;
    end; 

    if isempty(timePick) == 0 
        for iwave = 1:size(timePick, 2)
            if isempty(timePick(itrc, iwave)) == 0
                plot(timePick(itrc, iwave), recCoor(itrc), ...
                     timePickMarkerType, 'LineWidth', timePickLineWidth, ...
                     'MarkerSize', timePickMarkerSize, ...    %  'Color', timeCompColor);  
                     'MarkerEdgeColor', timePickColor, 'MarkerFaceColor', timePickColor);    
            end;  
        end;
    end;
    
end;    % of loop over 'itrc'
drawnow;

end    % of the function 