%% *plotSonic*
% Plot sonic log

%%
% *Input:*

% fig     - [scalar] number of a new figure
% inpData - [length(depth), :] array of tracks to be plotted
% depth   - [length(depth), 1] array of depths
% minX    - [scalar] - min limit of the X-axis
% maxX    - [scalar] - max limit of the X-axis
% minY    - [scalar] - min limit of the Y-axis
% maxY    - [scalar] - max limit of the Y-axis

%%
% *Output:* none

%%
% *Author:* Vladimir Grechka 2012

%%
function [] = plotSonic(fig, inpData, depth, minX, maxX, minY, maxY)
%% Settings
[~, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

lightGray  = 0.5*[1 1 1];   
colorArray = [[0 0 0]; lightGray; [0 0 1]; [1 0 0]];

%% Plot a log
figure(fig);
for icol = 1:size(inpData, 2)
    plot(inpData(:,icol), depth, '-', 'LineWidth', 0.5, 'Color', colorArray(icol,:));
end;

if isempty(minX) == 1;   minX = min(min(inpData));   end;
if isempty(maxX) == 1;   maxX = max(max(inpData));   end;
if isempty(minY) == 1;   minY = min(depth);          end;
if isempty(maxY) == 1;   maxY = max(depth);          end;
axis([minX, maxX, minY, maxY]);
drawnow;

end    % of the function