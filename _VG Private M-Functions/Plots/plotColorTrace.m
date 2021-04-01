%% *plotColorTrace*
% Plot a seismic trace with colored peaks and troughs  

%%
% *Input:*

% t           - [1, :] array of time samples; might be empty
% y           - [1, :] array of the trace values to be plotted
% lineColor   - [1, 3] color array for line to plot the trace
% lineWidth   - [scalar] width of line to plot the trace
%               (*) if lineWidth <= 0, no line is plotted 
% peakColor   - [1, 3] color array for plotting the peaks
% troughColor - [1, 3] color array for plotting the troughs
% alphaNum    - [scalar] number determining transparancy in Matlab function 'alpha'
%               (*) alphaNum < 1 is useful to display overlapping traces but time-consuming
%               (*) to speed up plotting, set alphaNum = [] or alphaNum = 1.0 

%%
% *Output:* none

%%
% *Author:* Vladimir Grechka 2012
%
% * function |plotColorTrace| with four or fewer input parameters plots wiggles

%%
function plotColorTrace(t, y, lineColor, lineWidth, peakColor, troughColor, alphaNum)
%% Settings
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

if isempty(y) == 1 
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Input trace array is empty \n \n');    
      error('>>> STOP');
end;

%% Plot colored trace
if isempty(peakColor) == 0  ||  isempty(troughColor) == 0   
    yBase = mean(y);
    
    % Complement the trace with its zero crossings to produce exact coloring of the peaks and 
    % troughs 
    t1 = t;    dt = t(2) - t(1);
    y1 = y - yBase;
    y2 = [y1(1), y1(1:end-1)];                          % shift the trace by one sample
    indX0 = find(y1.*y2 < 0);                           % find samples of the trace zero crossings
    t1(indX0) = t1(indX0-1) - dt*y1(indX0-1)./(y1(indX0) - y1(indX0-1));    % exact zero crossings
    
    count = 0;
    for i = 1:length(t)
        if t1(i) ~= t(i)
            count = count + 1;   tt(count) = t1(i);   yy(count) = yBase;
        end;
            count = count + 1;   tt(count) = t(i);    yy(count) = y(i);
    end;

    indPos = find(yy - yBase >  0);                     % find indexes of the positive and 
    indNeg = find(yy - yBase <= 0);                     % negative trace values
end;

if isempty(peakColor) == 0
    yPlot(indPos) = yy(indPos);                         % plot peaks
    yPlot(indNeg) = yBase;
    
    if isempty(yPlot) == 0
        area(tt, yPlot, yBase, 'FaceColor', peakColor, 'LineStyle', 'none');
        hold on;
    end;
end;

if isempty(troughColor) == 0
    yPlot(indNeg) = yy(indNeg);                         % plot troughs
    yPlot(indPos) = yBase;

    if isempty(yPlot) == 0
        area(tt, yPlot, yBase, 'FaceColor', troughColor, 'LineStyle', 'none');
        hold on;
    end;
end;

%% Apply transparency
if isempty(alphaNum) == 0
    alpha(alphaNum);                                
end;
                                                        
%% Plot a wiggle
if lineWidth > 0
    plot(t, y, 'LineStyle', '-', 'Color', lineColor, 'LineWidth', lineWidth);
end;

set(gca, 'YGrid', 'off');                           % turn off horizontal grid lines
set(gca, 'Layer', 'top');                           % bring vertical grid lines to the foreground 

end    % of the function
