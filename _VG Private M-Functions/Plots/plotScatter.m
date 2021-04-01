%% *plotScatter*
% Plot an array of data points in 3D

%%
% *Input:*

% fig          - [scalar] number of an existing figure
% array        - [3, :, :] array of the data points to be plotted
% inputUnits   - [string] input length units
% plotUnits    - [string] length units for plotting 
% symbolMarker - [char] marker used for plotting
% symbolSize   - [scalar] marker size  
% symbolColor  - [1, 3] color array for plotting

%%
% *Output:* none

%%
% *Author:* Vladimir Grechka 2012

%%
function plotScatter(fig, array, inputUnits, plotUnits, symbolMarker, symbolSize, symbolColor)
%% Settings  
[~, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

%% Plot data points to an existing figure
figure(fig); 
array = u2u(array, inputUnits, plotUnits, 0);           % change units 
for i = 1:size(array,2)
    plot3(squeeze(array(1,i,:)), squeeze(array(2,i,:)), squeeze(array(3,i,:)), ...
          symbolMarker, 'LineWidth', 0.5, ...
          'MarkerSize', symbolSize, 'MarkerEdgeColor', symbolColor, 'MarkerFaceColor', 'w'); 
end;
drawnow;

end    % of the function