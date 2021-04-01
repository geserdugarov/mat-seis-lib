%% *plotTraveltimes*
% Plot traveltimes

%%
% *Input:*

% fig             - [scalar] number of a new figure
% figPosition     - [4, 1] vector specified the figure position and size
% figTitle        - [string] figure 'title' 
% figXlabel       - [string] figure 'xlabel'
% figYlabel       - [string] figure 'ylabel'
%                   (*) The parameters above are the same as those for function 'setFigure' 
%                       if dim = 2
% time            - [noSou*noRec*noWave, 1] array of traveltimes
% tau0            - [noSou, 1] array of the origin times
% inputTimeUnits  - [char] input time units
% plotTimeUnits   - [char] time units for plotting 
% xRec            - [1, noRec] array of the receiver coordinates at which the times were
%                   measured or computed 
% inputLenthUnits - [string] input length units
% plotLengthUnits - [string] length units for plotting 
% markerSizeRec   - [scalar] marker size for receivers
% markerColorRec  - [1, 3] color array for receivers
% noSou           - [scalar] number of sources
% noRec           - [scalar] number of receivers
% noWave          - [scalar] number of wave codes
% markerSizeTime  - [scalar] marker size for traveltimes
% markerColorTime - [1, 3] color array for traveltimes
% recFlag         - [scalar] equal to 1 or 0 that controls whether or not to plot receivers  
% wellFlag        - [scalar] equal to 1 or 0 that controls whether or not to plot a well  
% eventName       - [cell] array of event names to be placed in a legend 
%                   (*) Empty 'eventName' causes replacing the event names with their 
%                       sequential numbers
% saveFigFlag     - [1, 3] array of 0 or 1 that controls saving a figure in 
%                   function saveFigure
% figNamePrefix   - [string] prefix of saved figures

%%
% *Output:*

% fignext         - [scalar] number of the next figure to be plotted

%%
% *Author:* Vladimir Grechka 2012 2013

%% 
% *Known issues:* 
%
% * 'legend' gives an error if noSymbol = 1

%%
function [fignext] = ...
      plotTraveltimes(fig, figPosition, figTitle, figXlabel, figYlabel, ...
                      time, tau0, inputTimeUnits,   plotTimeUnits, ...
                             xRec, inputLengthUnits, plotLengthUnits, ...
                      markerSizeRec,  markerColorRec,  noSou, noRec, noWave, ...
                      markerSizeTime, markerColorTime, recFlag, wellFlag, eventName, ...
                      saveFigFlag, figNamePrefix)
%% Settings
[~, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

waveCodeColor = setPlottingColors;
sourceSymbol  = ['s'; '<'; 'o'; '^'; 'd'; '>'; 'p'; 'v'; 'h'];
noColor  = size(waveCodeColor, 1);
noSymbol = size(sourceSymbol,1);

%% Change the units 
time = u2u(time, inputTimeUnits, plotTimeUnits, 0);
timePlot = reshape(time, noWave, noRec, noSou);
if isempty(tau0) == 0 
    tau0 = u2u(tau0, inputTimeUnits, plotTimeUnits, 0);
else
    tau0 = zeros(noSou, 1);
end;

xRec = u2u(xRec, inputLengthUnits, plotLengthUnits, 1);

%% Compute the times to be plotted
for isou=1:noSou
    timePlot(:,:,isou) = timePlot(:,:,isou) - tau0(isou);
end;

timePlotMin = min(min(min(timePlot)));
timePlotMax = max(max(max(timePlot)));

%% Plotting
figCount = 0;
for isou = 1:noSou   
    jsou = mod(isou, noSymbol);                                 % moveout number on a given plot
    
    if jsou == 0;   jsou = noSymbol;   end;
    if mod(isou, noSymbol) == 1
        if isou ~= 1
            fig = fig + 1;
        end;
        
        setFigure(fig, 2, figPosition, figTitle, figXlabel, figYlabel, [], [], []); 
        recPosition = timePlotMin*ones(size(xRec)) - 0.05*(timePlotMax - timePlotMin);
        if wellFlag == 1                                            % plot the well trajectory
            line('XData', recPosition, 'YData', xRec, 'Color', 0.75*[1 1 1], 'LineWidth', 3);
        end;
            
        if recFlag == 1                                             % plot the receiver locations
            plot(recPosition, xRec, '^', 'LineWidth', 0.5, 'MarkerSize', markerSizeRec, ...
                 'MarkerEdgeColor', 'w', 'MarkerFaceColor', markerColorRec); 
        end;

    end;

    % Plot the traveltimes
    for iwave = 1:noWave
        indexWave = mod(iwave, noColor);                % cycle the wave-type colors
        if indexWave == 0;   indexWave = noColor;   end;
    
        indexSou = mod(isou, noSymbol);                 % cycle the source symbols
        if indexSou == 0;   indexSou = noSymbol;   end;
        if isempty(markerColorTime) == 1;    
            markerColorTimePlot = waveCodeColor(indexWave,:);    
        else
            markerColorTimePlot = markerColorTime;        
        end; 

        markerSizeTimePlot = markerSizeTime;
        if sourceSymbol(indexSou) == 'p'  ||  sourceSymbol(indexSou) == 'h';
            markerSizeTimePlot = markerSizeTime + 2;
        end;

        plotHandle = plot(timePlot(iwave,:,isou), xRec, ...
            strcat('-', sourceSymbol(indexSou)), 'LineWidth', 0.5, ...
            'Color', waveCodeColor(indexWave,:), 'MarkerSize', markerSizeTimePlot, ...
            'MarkerFaceColor', markerColorTimePlot);

        if iwave == 1;
            % Create text for the legend  
            pltHandle(jsou) = plotHandle;    
            if isempty(eventName) == 1
                txt{jsou} = [' ', num2str(isou, '% 2.0f')]; 
            else
                tmp = eventName{isou};   
                txt{jsou} = [' ', num2str(isou, '% 2.0f'), ' - ', tmp]; 
            end;
        end;
    end;    % of loop over iwave

    if mod(isou, noSymbol) == 0  ||  isou == noSou
        lgnd = legend(pltHandle, txt, 'Location', 'SouthEastOutside');
        set(lgnd, 'Interpreter', 'none')
        clear txt; 
        axis('tight');
        drawnow;

        figCount = figCount + 1;
        if isempty(saveFigFlag) == 0
            saveFigure(fig, fullfile('Figures', ...
                  [figNamePrefix, '-traveltimes ', num2str(figCount, '%02.0f')]), saveFigFlag);
        end;
    end;
    
end;    % of loop over isou
fignext = fig;

end    % of the function