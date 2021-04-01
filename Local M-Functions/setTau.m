%% *setTau*
% Set the origin times 

%% 
% *NB:*
%
% * This script needs to be tailored to each specific data set

%%
% *Input:* 

% noSou - the number of sources

%%
% *Output:*

% tau  - [1, noSou] array of the origin times of perforation shots and microseismic events

%%
% *Author:* Vladimir Grechka 2012

%%
function [tau] = setTau(noSou)

% Settings  
% [thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
% narginTrue = nargin(thisFileName);   
% narginchk(narginTrue, narginTrue);

tau = zeros(1, noSou);
%tau = 0.1*(1 : noSou);

end  % of the function