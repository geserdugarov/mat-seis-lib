%% *linMap*
% Linear mapping of one array into another, typically to convert times to samples and vice versa

%% 
% *Input:*

% X  - input array to be linearly mapped that belongs to the scale [x0, xn]
% x0 - [scalar] value of the beginning of the scale for the input array
% dx - [scalar] increment of the scale for the input array
% y0 - [scalar] value of the beginning of the scale for the output array
% dy - [scalar] increment of the scale for the output array

%% 
% *Output:*

% Y  - output array, such that size(Y) = size(Y)

%%
% *Author:* Vladimir Grechka 2013

%%
function Y = linMap(X, x0, dx, y0, dy)
%% Settings 
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

if dx < eps
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Increment of the input array dx = %g cannot be that small \n', dx);
      error('>>> STOP');
end;

%% Linear mapping
Y = y0 + (dy/dx)*(X - x0);  
        
end  % of the function 