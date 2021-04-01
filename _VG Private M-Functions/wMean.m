%% *wMean*
% Weighted mean and std of an array computed along its first dimension 

%%
% *Input:*

% arrIn   - [:, :] input array  
% stdIn   - [:, :] array of standard deviations of elements of 'arrIn'

%%
% *Output:*

% meanOut - [1, :] weighted mean of 'arrIn'
% stdOut  - [1, :] weighted std of 'arrIn'

%%
% *Author:* Vladimir Grechka 2012

%%
function [meanOut, stdOut] = wMean(arrIn, stdIn)
%% Settings and defaults
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

stdInTol = stdIn;
tol = 1.e-12;    % tolerance to avoid division by zero

meanOut = NaN(1, size(arrIn,2));    
stdOut  = NaN(1, size(arrIn,2));    % initialize

%% Calculate meanOut and stdOut (see <http://en.wikipedia.org/wiki/Weighted_mean>)
for i2 = 1:size(arrIn,2) 
    % Find ~NaNs in both 'arrIn' and 'stdIn'
    notNaN = find(~isnan(arrIn(:,i2)) & ~isnan(stdIn(:,i2)));  
    
    if isempty(notNaN) == 0
        % There are data points for which both 'arrIn' and 'stdIn' are not NaNs
        stdInTol(stdIn(:,i2) < tol, i2) = tol;    % perturb near 0 std's to avoid division by 0

        smean = sum(arrIn(notNaN,i2)./stdInTol(notNaN,i2));     % Data points for which
        sstd  = sum(1./stdInTol(notNaN,i2));                    % stdIn = NaN are excluded from
        meanOut(1,i2) = smean/sstd;                             % the calculation of 'meanOut'
        
        if sstd == 1
            stdOut(1,i2) = sstd;     % avoids division by zero
            fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
            fprintf('>>> WARNING: Defaulting intermediate ''stdOut'' to %g to avoid singularity \n', ...
                sstd);
        else
            stdOut(1,i2) = ...
                sqrt(sum((arrIn(notNaN,i2) - meanOut(1,i2)).^2./stdInTol(notNaN,i2))/(sstd - 1));   
                % This formulation ensures that, for example,
                % wMean([2,2,2,4,5,5]', [1,1,1,1,1,1]') = wMean([2,4,5]', [1/3,1,1/2]')
        end;
        % Make sure that 'stdOut' is not smaller than the minimum of 'stdIn'
        stdOut(1,i2) = max([stdOut(1,i2), min(stdIn(:,i2))]);
        
    else
        % stdIn = NaN for all data -- calculate the unweighted mean and std
        notNaNarr = find(~isnan(arrIn(:,i2)));                      % ~NaN data but NaN std's
        if isempty(notNaNarr) == 0
            meanOut(1,i2) = mean(arrIn(notNaNarr,i2));              % replace 'wMean' with mean 
            stdOut (1,i2) =  std(arrIn(notNaNarr,i2));              % and 'wStd' with std
        end;
    end;
end;

end    % of the function