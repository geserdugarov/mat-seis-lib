%% *sortTraces*
% Sort the components of 3C seismic data to [East, North, Down]

%%
% *Input:*

% traceData  - [:, 3] array containing the input 3C trace
% headerData - [structure] of trace headers produced by function 'readSEGYSEG2'
% polarity   - [1, 3] array containing +1 or -1 to indicate the polarity flip of 'traceData'
% compOrder  - [1, 3] array that associates the components of (possibly polarity-flipped) array
%              'traceData' with the left-hand coordinate frame [East, North, Down]
%              (*) For example, sorting the components oriented as [Up, West, North] leads to 
%                  compOrder = [3, 1, 2] and polarity = [-1, -1, 1]

%%
% *Output:*

% traceOut   - [:, 3] array 'traceData' whose components are sorted to [East, North, Down]      

%%
% *Author:* Vladimir Grechka 2012

%%
function [traceOut] = sortTraces(traceData, headerData, polarity, compOrder)
%% Settings and defaults
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

traceOut = NaN(size(traceData));    tmp1 = NaN(size(traceData, 1), 3);    tmp2 = tmp1;

if isempty(traceData) == 1
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Empty trace array \n');    
    fprintf('>>> PAUSE -- Continue? \n');    
    return;
end;

%% Change the polarity of components and resort them
if isequal(abs(polarity), [1, 1, 1]) ~= 1
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Improper array polarity = [%g, %g, %g] \n', polarity);
    fprintf('>>> Array ''polarity'' should contain +1 and -1 \n \n');
      error('>>> STOP');    
elseif isequal(sort(compOrder), [1, 2, 3]) ~= 1
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Improper array compOrder = [%g, %g, %g] \n', compOrder);
    fprintf('>>> Array ''compOrder'' should contain a permutation of 1, 2, 3 \n \n');
      error('>>> STOP');    
else
    for irec = 1:headerData.noRec
        tmp1 = traceData(:, 3*irec-2:3*irec).*repmat(polarity, size(traceData, 1), 1);
        tmp2(:, compOrder) = tmp1;
        traceOut(:, 3*irec-2:3*irec) = tmp2;
    end;
end;    

end    % of the function