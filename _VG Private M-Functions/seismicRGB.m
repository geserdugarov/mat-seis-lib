%% *seismicRGB*
% Create seismic colormap, , ranging from dark blue via white to dark red

%%
% *Input:* none

%%
% *Output:* none

%%
% *Author:* Nico Sneeuw, Munich, 31/08/94

% Original script, called 'seismic' can be downloaded from
% 'http://www.mathworks.com/matlabcentral/fileexchange/
% 30585-large-data-in-matlab--a-seismic-data-processing-case-study/content/migration/seismic.m'

% function rgb = seismic(n)
% seismic(n) creates a colormap, ranging from dark blue via white to dark red.
% Nico Sneeuw
% Munich, 31/08/94

%%
function rgb = seismicRGB
%% Settings 
[~, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

%if nargin == 0, 
    n = size(get(gcf,'colormap'),1); 
%end

%% Create colormap
m = ceil(n/3);
top = ones(m,1);
bot = zeros(m,1);
up = (0:m-1)'/m;
down = flipud(up);

r = [bot; up; 1; top; down];
g = [bot; up; 1; down; bot];
b = [up; top; 1; down; bot];
rgb = [r g b];

% rgb-map has size 4m+1 now. The central part will be extracted.
xlarge = 4*m+1-n;
xblue = round(xlarge/2);
xred = xlarge - xblue;
rgb([1:xblue 4*m-xred+2:4*m+1],:) = [];

end    % of the function
