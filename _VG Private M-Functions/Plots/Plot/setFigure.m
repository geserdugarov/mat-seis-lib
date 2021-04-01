%% *setFigure*
% Set up a new figure

%% 
% *Input:*

% fig         - [scalar] number of a new figure 
% dim         - [scalar] dimension of the figure; 'dim' should be equal to 2 or 3  
% figPosition - [1, 4] vector specifying the figure position and size in the form 
%                      [left, bottom, width, height]
% figTitle    - [string] figure 'title' 
% figXlabel   - [string] figure 'xlabel'
% figYlabel   - [string] figure 'ylabel'
% figZlabel   - [string] figure 'zlabel'
% figFontSize - [scalar] font size
% titleInterp - [string] equal to either 'latex' (default), 'tex', or 'none' to turn on and off
%               LaTeX or TeX interpretation of the title string
% mn          - [1, 2] vector of the 'm-n' parameters of Matlab function 'subplot'
% p           - [scalar or vector] of the 'p' parameter of Matlab function 'subplot'
%               (*) Empty or absent parameters 'mn' and 'p' invoke Matlab 'figure' without subplots

%% 
% *Output:* none

%% 
% *Author:* Vladimir Grechka 2012 - 2014

%%
function setFigure(fig, dim, figPosition, figTitle, figXlabel, figYlabel, figZlabel, ...
                   figFontSize, titleInterp, mn, p)
%% Settings and defaults
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
%narginTrue = nargin(thisFileName);   
%narginchk(narginTrue, narginTrue);

if ~(dim == 2  ||  dim == 3)  
    fprintf('>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Incorrect dimension of a new figure dim = %g \n', dim);
    fprintf('>>> The correct values are dim = 2 or dim = 3 \n \n');  
      error('>>> STOP');
end;

if isempty(figPosition) == 1
    scrsz = get(0, 'ScreenSize');
    corner = 1/6;
    figPosition = [corner*scrsz(3), corner*scrsz(4), ...
           (1 - 2*corner)*scrsz(3), (1 - 2*corner)*scrsz(4)];
end;

if isempty(figFontSize) == 1;   figFontSize = 16;   end;
 
if (isequal(titleInterp, 'latex') + isequal(titleInterp, 'tex') + ...
    isequal(titleInterp, 'none')) == 0   
    titleInterp = 'tex';   
end;

%% Initialize a figure or a subplot 
figure(fig);  
if nargin == 9
    % Initialize a figure
    set(fig, 'Color', [1 1 1]);     % set color for the figure window
    set(fig, 'Position', figPosition); 
elseif isempty(mn) == 0  &&  isempty(p) == 0
    % Initialize a subplot
    subplot(mn(1), mn(2), p);
end;

%% Define properties of the figure or the subplot
set(gca, 'FontSize', figFontSize, 'FontName', 'Times New Roman'); 
grid on;  hold on; 

if dim == 2
    box on;                     % plot border 
    set(gca, 'Ydir', 'reverse'); 
elseif dim == 3
    box off;
    set(gca, 'Zdir', 'reverse'); 
    zlabel(figZlabel);   
    view(3);   
end;

xlabel(figXlabel);   ylabel(figYlabel); 
title(figTitle, 'FontSize', (figFontSize + 2), 'FontName', 'Times New Roman', ...
    'Interpreter', titleInterp);
drawnow;

end   % of the function 
