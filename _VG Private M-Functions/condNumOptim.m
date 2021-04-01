%% *condNumOptim*
% Objective function for obtaining the scaling that minimizes the condition number

%%
% *Input:*

% param     - [1, 3] vector of scales of the units of time, distance, and angle -- in that order
% Mat       - rectangular matrix whose condition number is to be minimized 
% expUnits  - [3, size(Mat,2)] matrix whose rows contain the exponentials of time, distance, and  
%             angle so the units of the i-th column are 
%             T^expUnits(1,i)*L^expUnits(2,i)*A^expUnits(3,i)
% printFlag - [scalar] flag, equal to 0 or 1, to indicate verbose print 

%%
% *Output:*

% F        - [scalar] log10 of the condition number

%%
% *Author:* Vladimir Grechka 2012 - 2014

%%
function [F] = condNumOptim(param, Mat, expUnits, printFlag)

global scaleMat scaleCond numCond

%% Settings and defaults
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

%% Compute the condition number
scaleMat = Mat;    scaleCond = ones(1,size(Mat,2));
for i=1:size(Mat,2)                                                 
    scaleCond(1,i) = 1/prod(param(:).^expUnits(:,i), 1);  % apply the scaling factors
    scaleMat(:,i)  = Mat(:,i)*scaleCond(1,i);             % multiply columns of Mat by 'scaleCond'
end;  
numCond = cond(scaleMat);

%% The objective function
F = log10(numCond);                                 

% [log10(min(scaleCond)), log10(max(scaleCond)), F]

%% Output and warning messages
if printFlag == 1
    fprintf('>>> Log10 of condition number = %g \n', F);
end;

if isnan(numCond) == 1  ||  isinf(numCond) == 1
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> WARNING: Condition number = %g -- PAUSE \n', numCond);    pause
end;

end    % of the function