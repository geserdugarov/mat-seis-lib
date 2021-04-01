%% *setPlottingColors*
% Set colors to mark the populations of ray trajectories and moveouts 

%%
% *Input:* none

%%
% *Output:*

% rayColors - [:, 3] array of colors  

%%
% *Author:* Vladimir Grechka 2013

%%
function [rayColors] = setPlottingColors
%% Settings 
% [~, thisFileName, ~] = fileparts(mfilename('fullpath'));
% narginTrue = nargin(thisFileName);   
% narginchk(narginTrue, narginTrue);

rayColors  = [ ...
                  [0 0 0]; ...              % black 
                  [1 0 0]; ...              % red 
                  [0 0 1]; ...              % blue
                  [0 1 1]; ...              % cyan
                  [0 1 0]; ...              % green
                  [1 0 1]; ...              % magenta
              0.5*[1 1 1]; ...              % gray 
                  [1 1 0]  ...              % yellow
              ];

end    % of the function