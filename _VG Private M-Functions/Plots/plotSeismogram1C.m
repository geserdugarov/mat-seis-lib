%% *plotSeismogram1C*
% Plot single-component wiggle seismogram

%%
% *Input:*

% fig                - [scalar] number of a new figure
% time               - [1, noTime] array of time samples
%                      (*) if array time is empty, the traces are ploted in samples
% recCoor            - [noRec, 1] - array of receiver coordinates 
%                      (*) if array recCoor is empty, the seismogram is labeled in trace numbers
% traceDataT         - [noTime, noRec] matrix containing traces
% flagNorm           - [scalar] normalization flag
%                      flagNorm = 0 - no normalization (default)
%                      flagNorm = 1 - normalize to the maximim of the entire seismogram
%                      flagNorm = 2 - normalize each trace individually
% ampFactor          - [scalar] amplification factor for plotting wiggle traces (style = 'wiggle')
%                      (*) if variable ampFactor is empty, ampFactor is defaulted to 1
% clipLevel          - [scalar] from 0 to 1 (default) determining the clip level of seismic
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
% traceLineColor     - [1, 3] color array for lines to plot the traces (see 'plotColorTrace')
% traceLineWidth     - [scalar] width of lines to plot the traces (see 'plotColorTrace')
%                      (*) if traceLineWidth <= 0, no line is plotted 
% peakColor          - [1, 3] color array for plotting the peaks (see 'plotColorTrace')
% troughColor        - [1, 3] color array for plotting the troughs (see 'plotColorTrace')
% alphaNum           - [scalar] number determining transparancy (see 'plotColorTrace')
%                      (*) alphaNum < 1 is useful to display overlapping traces but might be 
%                          time consuming
%                      (*) to speed up plotting, set alphaNum = [] or alphaNum = 1.0 
% orientation        - [string] indicating the orientation of the time axis:
%                      'vertical' or 'horizontal'
% style              - [string] indicating the style of display: 'wiggle' or 'image'

%%
% *Comments:*
%
% * If the parameter orientation = 'horizontal', traces cannot be colored and the parameters 
%   'peakColor', 'troughColor', and 'alphaNum' become irrelavant. This limitation stems from the 
%   features of the Matlab function 'area' that can color areas up or down but not left or right
%
% * Also parameters 'ampFactor', 'traceLineColor', 'traceLineWidth', 'peakColor', 'troughColor', 
%   and 'alphaNum', controlling the wiggle display, are irrelevant when the parameter 'style' is
%   set to 'style' = 'image'

%%
% *Output:* none

%%
% *Author:* Vladimir Grechka 2012 - 2014

%%
function plotSeismogram1C(fig, time, recCoor, traceDataT, ...
    flagNorm, ampFactor, clipLevel, polarity, tmin, tmax, ...
    timePick, timePickColor, timePickMarkerType, timePickLineWidth, timePickMarkerSize, ...
    timeComp, timeCompColor, timeCompMarkerType, timeCompLineWidth, timeCompMarkerSize, ...
    traceLineColor, traceLineWidth, peakColor, troughColor, alphaNum, orientation, style)
%% Settings 
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

traceData = traceDataT';
[noRec, noTime] = size(traceData);
wiggleDisplayLevel = 1.0;

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

if isempty(flagNorm) == 1  ||  sum(flagNorm == [0, 1, 2]) == 0
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

if isempty(clipLevel) == 1   
    clipLevel = 1;                                  % default for variable 'clipLevel'
elseif clipLevel < 0  ||  clipLevel > 1
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf(['>>> WARNING: Variable ''clipLevel'' [= %g] is defaulted to 1 because lies ', ...
             'outside the expected interval [0, 1] \n'], clipLevel);    
    display('>>> PAUSE');  beep on;  beep;  beep off;  pause;  
    clipLevel = 1; 
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

if isempty(orientation) == 1
    orientation = 'horizontal';                     % default for variable 'orientation'
end;

if isempty(style) == 1
    style = 'wiggle';                               % default for variable 'style'
end;

%% Set up plotting
figure(fig);
grid on;   
if strcmp(orientation, 'horizontal') == 1 
    set(gca, 'YGrid', 'off');   
elseif strcmp(orientation, 'vertical') == 1 
    set(gca, 'XGrid', 'off');   
else
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Improper value of input parameter ''orientation'' = %s \n', orientation);    
    fprintf('>>> The correct values are ''vertical'' or ''horizontal'' \n \n');    
      error('>>> STOP');
end;
hold on; 

if numel(recCoor) > 1 
    dCoor = abs(recCoor(2) - recCoor(1));
else
    dCoor = 1;
end;
if strcmp(orientation, 'horizontal') == 1 
    xlim([tmin, tmax]);
    ylim([min(recCoor) - 0.999*dCoor, max(recCoor) + 0.999*dCoor]);
elseif strcmp(orientation, 'vertical') == 1 
    ylim([tmin, tmax]);
    xlim([min(recCoor) - 0.999*dCoor, max(recCoor) + 0.999*dCoor]);
end;

if polarity == 1
    if strcmp(orientation, 'horizontal') == 1 
        set(gca, 'YDir', 'normal'); 
    elseif strcmp(orientation, 'vertical') == 1 
        set(gca, 'XDir', 'normal'); 
    end;
elseif polarity == -1;
    if strcmp(orientation, 'horizontal') == 1 
        set(gca, 'YDir', 'reverse');      
    elseif strcmp(orientation, 'vertical') == 1 
        set(gca, 'XDir', 'reverse');      
    end;
else
    set(gca, 'YDir', 'normal');    
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> WARNING: Incorrect value of polarity = %g \n', polarity);    
    fprintf('>>> Variable ''polarity'' is defaulted to 1 \n \n');    
end;

%% Normalize the traces
data = reshape(traceData, noRec*noTime, 1);
maxData = max(abs(data));
if flagNorm == 1  &&  maxData > 0
    data = data/maxData;                            % normalize the entire seismogram
    traceData = reshape(data, noRec, noTime);
end;

%% Clip the traces
if strcmp(style, 'image') == 1 
    if flagNorm < 2
        data(data >  clipLevel*maxData) =  clipLevel*maxData;    
        data(data < -clipLevel*maxData) = -clipLevel*maxData;
        traceData = reshape(data, noRec, noTime);
    end;
end;

%% Plot the seismogram
for itrc = 1:noRec
    trc = traceData(itrc,:); 
    
    if flagNorm == 2
        maxData = max(abs(trc));
        if maxData > 0
            trc = trc/maxData;                      % normalize each trace individually
        end;

        if strcmp(style, 'image') == 1              % apply clipping
            trc(trc >  clipLevel*maxData) =  clipLevel*maxData;    
            trc(trc < -clipLevel*maxData) = -clipLevel*maxData;
            traceData(itrc,:) = trc;                % save the traces for image display
        end;
    end;
    
    if strcmp(style, 'wiggle') == 1                 % plot the traces as wiggles
        trc = ampFactor*trc;                        % apply the amplitude factor
        trc(trc >  wiggleDisplayLevel) =  wiggleDisplayLevel;    % clip the traces for plotting
        trc(trc < -wiggleDisplayLevel) = -wiggleDisplayLevel;

        if strcmp(orientation, 'horizontal') == 1 
%        plot(time, recCoor(itrc) + dCoor*trc, ...
%            'LineStyle', '-', 'Color', traceLineColor, 'LineWidth', traceLineWidth);
            plotColorTrace(time, recCoor(itrc) + dCoor*trc, ...
                traceLineColor, traceLineWidth, peakColor, troughColor, alphaNum);
        elseif strcmp(orientation, 'vertical') == 1 
            plot(recCoor(itrc) + dCoor*trc, time, ...
                'LineStyle', '-', 'Color', traceLineColor, 'LineWidth', traceLineWidth);
        end;
    end;
end;

if strcmp(style, 'image') == 1                      % plot the traces as image
    if strcmp(orientation, 'horizontal') == 1 
%        pcolor(time, recCoor, traceData);
        imagesc(time, recCoor, traceData);
    elseif strcmp(orientation, 'vertical') == 1 
        imagesc(recCoor, time, traceData');
    end;
%    colormap(gray);
%    colormap(seismic(3));                           % Mauricio Sacchi's seismic colormap
    colormap(seismicRGB);                           % Matlab seismic colormap
    set(gca, 'Layer', 'top');                       % bring up the grid lines to the top
end;

%% Plot the computed and picked times
for itrc = 1:noRec
    if isempty(timeComp) == 0
        for iwave = 1:size(timeComp, 2)
            if isempty(timeComp(itrc, iwave)) == 0
                if strcmp(orientation, 'horizontal') == 1 
                    plot(timeComp(itrc, iwave), recCoor(itrc), ...
                         timeCompMarkerType, 'LineWidth', timeCompLineWidth, ...
                         'MarkerSize', timeCompMarkerSize, ... % 'Color', timeCompColor);  
                         'MarkerEdgeColor', timeCompColor, 'MarkerFaceColor', timeCompColor);
                elseif strcmp(orientation, 'vertical') == 1 
                    plot(recCoor(itrc), timeComp(itrc, iwave), ...
                         timeCompMarkerType, 'LineWidth', timeCompLineWidth, ...
                         'MarkerSize', timeCompMarkerSize, ... % 'Color', timeCompColor);  
                         'MarkerEdgeColor', timeCompColor, 'MarkerFaceColor', timeCompColor);
                end;
            end;       
        end;
    end; 

    if isempty(timePick) == 0 
        for iwave = 1:size(timePick, 2)
            if isempty(timePick(itrc, iwave)) == 0
                if strcmp(orientation, 'horizontal') == 1 
                    plot(timePick(itrc, iwave), recCoor(itrc), ...
                         timePickMarkerType, 'LineWidth', timePickLineWidth, ...
                         'MarkerSize', timePickMarkerSize, ...      %  'Color', timeCompColor);  
                         'MarkerEdgeColor', timePickColor, 'MarkerFaceColor', timePickColor);    
                elseif strcmp(orientation, 'vertical') == 1 
                    plot(recCoor(itrc), timePick(itrc, iwave), ...
                         timePickMarkerType, 'LineWidth', timePickLineWidth, ...
                         'MarkerSize', timePickMarkerSize, ...      %  'Color', timeCompColor);  
                         'MarkerEdgeColor', timePickColor, 'MarkerFaceColor', timePickColor);    
                end;
            end;  
        end;
    end;
    
end;    % of loop over 'itrc'
drawnow;

end    % of the function

%%
