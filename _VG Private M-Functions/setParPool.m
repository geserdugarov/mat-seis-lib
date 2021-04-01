%% *setParPool*
% Set up parallel computations

%% 
% *Input:* none

%% 
% *Output:* 

% noWorkers - [scalar] number of cored for parallel computations


%% 
% *Author:* Vladimir Grechka 2014

%%
function [noWorkers] = setParPool
%% Settings and defaults
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

%% Determine the Matlab version and start parallel computations
verMatlab  = ver('symbolic');
yearMatlab = str2num(verMatlab.Release(3:end-2));
if yearMatlab < 2014
    % Older Matlab versions
    mps = matlabpool('size');  
    if mps == 0
        matlabpool open local
    end;
    noWorkers = matlabpool('size');  
else
    % The R2014a Matlab version and newer
    pool = gcp('nocreate');
    if isempty(pool) == 1
        parpool;
    end;
    pool = gcp;
    noWorkers = pool.NumWorkers;
end;

fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
fprintf('>>> %g workers for parallel computing found \n \n', noWorkers);

end   % of the function 
