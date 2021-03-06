%% *getInterface*
% Get interfaces for the inversion

%%
% *Input:*

% unknInd   - [5, noInt] index array generated by function 'setUnknowns'
% param     - [1, :] array of the intrface-related unknowns

%%
% *Output:*

% interface - [5, noInt] array obtained from function 'setInterface' and subsequently modified
% noInt     - the number of interfaces

%%
% *Author:* Vladimir Grechka 2012 

%%
function [interface, noInt] = getInterface(unknInd, param)
%% Settings 
[~, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

[interface, noInt] = setInterface;

%% Modify elements of array 'interface'
[irow, icol] = find(unknInd == 0);
for i = 1:length(irow)
    interface(irow(i), icol(i)) = param(i);                 
end;

end    % of the function