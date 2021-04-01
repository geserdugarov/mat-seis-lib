%% *plotBars*
% Plot error bars 

%%
% *Input:*

% fig           - [scalar] number of an existing figure
% yArray        - [:, :] array of quantities used to generate the bars
% yConfInt      - [scalar] for the size of confidence interval
%                 Typicaly, yConfInt = 1 for a 68% confidence interval and yConfInt = 2 for 95%
% yArrayTrue    - [size(yConfInt,1), 1] array of the true values of yArray, if available
% yMin          - [scalar] minimum value for plotting
% yMax          - [scalar] maximum value for plotting
% yTick         - [1, :] array of Y-ticks 
% inputUnits    - [string] input units
% plotUnits     - [string] units for plotting 
% xCapt         - [cell] array for captions of the x-axis
% xCaptArray    - [1, :] array of the numbers of cells of 'xCapt' to be displayed
% xCaptFontSize - [scalar] caption font size
% markSymb      - [char] marker for plotting the mean values
% sizeSymb      - [scalar] marker size for the mean values
% colorSymb     - [1, 3] color array for plotting the bars
% markSymbTrue  - [char] marker for plotting the true values
% sizeSymbTrue  - [scalar] marker size for the true values 
% widthSymbTrue - [scalar] line width for plotting the true values
% colorSymbTrue - [1, 3] color array for plotting the true values

%%
% *Output:* none

%%
% *Author:* Vladimir Grechka 2012 - 2014

%%
function plotBars(fig, yArray, yConfInt, yArrayTrue, ...
                       yMin, yMax, yTick, inputUnits, plotUnits, ...
                       xCapt, xCaptArray, xCaptFontSize, ...
                       markSymb, sizeSymb, colorSymb, ...
                       markSymbTrue, sizeSymbTrue, widthSymbTrue, colorSymbTrue)
%% Settings and checks
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

if ismatrix(yArray) == 0
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Input array ''yArray'' is not a matrix \n \n'); 
      error('>>> STOP');
end;

if isempty(inputUnits) == 0  &&  isempty(plotUnits) == 0
   yArray     = u2u(yArray, inputUnits, plotUnits);         % change the units
   yArrayTrue = u2u(yArrayTrue, inputUnits, plotUnits);
   yMin       = u2u(yMin, inputUnits, plotUnits);
   yMax       = u2u(yMax, inputUnits, plotUnits);
   yTick      = u2u(yTick, inputUnits, plotUnits);
end;

if isempty(yTick) == 0;    set(gca, 'YTick', yTick);    end;

xArray = (1:size(yArray,1));
xlim([0.5, xArray(end)+0.5]);

%% Error bar plot
figure(fig);
if size(yArray, 2) == 1
    yMean = yArray;
    yErr  = zeros(size(yMean));
    plot(xArray, yArray, markSymb, 'LineWidth', 0.5, ...
        'MarkerEdgeColor', colorSymb, 'MarkerFaceColor', colorSymb);
else
    yMean = mean(yArray, 2);
    yErr  =  std(yArray, 0, 2);
    errorbar(xArray, yMean, yConfInt*yErr,  strcat(markSymb, colorSymb), ...
        'LineWidth', 2, 'MarkerSize', sizeSymb, 'MarkerEdgeColor', colorSymb); 
end;

if isempty(yArrayTrue) == 0
    plot(xArray, yArrayTrue, markSymbTrue, 'LineWidth', widthSymbTrue, ...
        'MarkerSize', sizeSymbTrue, 'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', colorSymbTrue); 
end;

% Size the plot
yBeg = min([yMin, yMean' - yConfInt*yErr', yArrayTrue']);
yEnd = max([yMax, yMean' + yConfInt*yErr', yArrayTrue']);
dy = 0.05*(yEnd - yBeg);
if dy > eps
    ylim([yBeg - dy, yEnd + dy]);
end;

%% Label the parameters
set(gca, 'XTick', xArray);
set(gca, 'XTickLabel', []);
for i = 1:length(xArray)
    text(i, yBeg - 1.5*dy, ...
        xCapt{xCaptArray(i)}, 'Fontsize', xCaptFontSize, 'Fontname', 'Times', ... 
        'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Top');
    % strcat(xCapt{1,xCaptArray(i)}(1:10),'}'), 'Fontsize', xCaptFontSize, 'Fontname', 'Times', ...
    % -- caption without the layer numbers, if desirable
end;
drawnow;

end    % of the function