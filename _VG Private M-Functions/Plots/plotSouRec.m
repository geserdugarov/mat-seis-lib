%% *plotSouRec*
% Plot the positions of sources and receivers

%%
% *Input:*

% fig             - [scalar] number of an existing figure
% xSou            - [3, :] array of the source coordinates
% xRec            - [3, :] array of the receiver coordinates
% inputUnits      - [string] input length units
% plotUnits       - [string] length units for plotting 
% indexPerf       - [1, size(xSou, 2)] array of the indexes of perforation shots specified in
%                   function 'setSouRec'
% markerSizePerf  - [scalar] marker size for perforation shots 
% markerColorPerf - [1, 3] color array for perforation shots
% indexEvnt       - [1, size(xSou, 2)] array of the indexes of microseismic events specified in
%                   function 'setSouRec'
% markerSizeEvnt  - [scalar] marker size for microseismic events
% markerColorEvnt - [1, 3] color array for microseismic events
% markerSizeRec   - [scalar] marker size for receivers
% markerColorRec  - [1, 3] color array for receivers
% wellFlag        - [scalar] equal to 1 or 0 that controls whether or not to plot a well  
%                   trajectory
% eventFlag       - [scalar] equal to 1 or 0 that controls whether or not to display the event
%                   numbers  

%%
% *Output:* none

%%
% *Author:* Vladimir Grechka 2012 2013

%%
function plotSouRec(fig, xSou, xRec, inputUnits, plotUnits, ...
                    indexPerf, markerSizePerf, markerColorPerf, ...
                    indexEvnt, markerSizeEvnt, markerColorEvnt, ...
                               markerSizeRec,  markerColorRec, wellFlag, eventFlag)
%% Settings

% [~, thisFileName, ~] = fileparts(mfilename('fullpath'));
% narginTrue = nargin(thisFileName);   
% narginchk(narginTrue, narginTrue);

xSou = u2u(xSou, inputUnits, plotUnits, 0);               % change units 
xRec = u2u(xRec, inputUnits, plotUnits, 0);

%% Plot positions of the sources and receivers
figure(fig);

for isou = 1 : size(xSou, 2)
    if indexPerf(isou) == 1
        % Plot perforation-shot locationss
        plot3(xSou(1,isou), xSou(2,isou), xSou(3,isou), 'p', 'LineWidth', 0.5, ...
              'MarkerSize', markerSizePerf, 'MarkerEdgeColor', 'w', ...
              'MarkerFaceColor', markerColorPerf);
    elseif indexEvnt(isou) == 1
        % Plot event locations
        plot3(xSou(1,isou), xSou(2,isou), xSou(3,isou), 'o', 'LineWidth', 0.5, ...
              'MarkerSize', markerSizeEvnt, 'MarkerEdgeColor', 'w', ...
              'MarkerFaceColor', markerColorEvnt);
%              'MarkerSize', markerSizeEvnt, 'MarkerEdgeColor', markerColorEvnt, ...
    end;
    
    if nargin > 13  &&  eventFlag == 1
        % Label perforation shots and events
        text(xSou(1,isou), xSou(2,isou), xSou(3,isou), ...
             strcat(' _{ }', num2str(isou)), 'Fontsize', 10, ...
             'Fontname', 'Times', 'Fontweight', 'Bold', ...
             'HorizontalAlignment', 'Left', 'VerticalAlignment', 'Middle');

    end;
end;

%% Plot the well trajectory
if nargin > 12  &&  wellFlag == 1  
    whiteColor = [1 1 1];   lightGrayColor = 0.75*whiteColor;
    if isempty(xRec) == 0
        line('XData', xRec(1,:), 'YData', xRec(2,:), 'ZData', xRec(3,:), ...
             'Color', lightGrayColor, 'LineWidth', 4);
    end;
end;

if isempty(xRec) == 0
    plot3(xRec(1,:), xRec(2,:), xRec(3,:), '^', 'LineWidth', 0.5, ...
        'MarkerSize', markerSizeRec, 'MarkerEdgeColor', 'w', 'MarkerFaceColor', markerColorRec);
%    plot3(xRec(1,:), xRec(2,:), xRec(3,:), 'o', 'LineWidth', 0.5, ...
%        'MarkerSize', markerSizeRec, 'MarkerEdgeColor', markerColorRec, 'MarkerFaceColor', markerColorRec);
end;

axis('image');   axis('tight');   
view(3);   
drawnow;

end    % of the function