%% *figureGrid*
% Create a grid for placing figures

%%
% *Input:*

% dimGrid  - [1, 2] array defining the number of rows and colums of the grid
% lrPixels - [scalar] offsets from the left and right edges of the screen (in pixels)

%%
% *Output:*

% figGrid  - [prod(dimGrid), 4] matrix whose rows are the figure [left, bottom, width, height]
%            as used in set(fig, 'Position', [left, bottom, width, height])

%%
% *Author:* Vladimir Grechka 2013

%%
function [figGrid] = figureGrid(dimGrid, lrPixels)
%% Settings and defaults
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
%narginTrue = nargin(thisFileName);   
%narginchk(narginTrue, narginTrue);

figGrid = NaN(prod(dimGrid), 4);
noRow = dimGrid(1);  noCol = dimGrid(2);

if noRow < 1 || noCol < 1
    fprintf('>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Incorrect dimensions of the figure grid = [%g, %g] \n', dimGrid);
    fprintf('>>> They should be integers greater or equal than 1 \n \n');  
      error('>>> STOP');
end;
    
% Defaults for the geometric parameters
bottomTaskbar = 50;  figBorder = 15;  figToolbar = 75;  figSpace = 15;  
if nargin < 2 || isempty(lrPixels) == 1;  lrPixels = 10;  end;  

%% Define 'figGrid'
scrsz = get(0, 'ScreenSize');    
figX = floor( (scrsz(3) - 2*lrPixels)/noCol - figBorder - figSpace ); 
figY = floor( (scrsz(4) - bottomTaskbar)/noRow - figBorder - figSpace - figToolbar  );

ifig = 0;
for icol = 1:noCol
    for irow = noRow:-1:1
        ifig = ifig + 1;
        figGrid(ifig, :) = ...
            [(icol - 1)*(figX + figBorder + figSpace) + lrPixels + 1, ...
             (irow - 1)*(figY + figBorder + figSpace + figToolbar) + bottomTaskbar + 1, ...
             figX, figY];
%        figure(ifig);  set(ifig, 'Position', figGrid(ifig,:));  title(num2str(ifig));  
%        drawnow;  fprintf('>>> PAUSE \n');  pause;   % uncomment to see instant action
    end;
end;

end  % of the function

%%
