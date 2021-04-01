%% *wAzim*
% Weighted |mean| and |std| of an array of azimuths computed along its first dimension 

%%
% *Input:*

% azmIn      - [:, :] array of the azimuths  
% stdIn      - [:, :] array of standard deviations of azimuths 
% angleUnits - [char] angle units of input variables 'azmIn' and 'azmIn' and output variables 
%              'meanOut', 'stdOut'
% tol        - [scalar] tolerance for issuing a warning message

%%
% *Output:*

% meanOut - [:, 1] weighted mean azimuth
% stdOut  - [:, 1] weighted std of azmIn

%%
% *Author:* Vladimir Grechka 2012

%%
function [meanOut, stdOut] = wAzim(azmIn, stdIn, angleUnits, tol)
%% Settings and defaults
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

if isempty(tol) == 1
    tol = 1.e-2;   
%    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
%    fprintf('>>> WARNING: Empty variable ''tol'' is defaulted to %g \n', tol);
end;

nData = size(azmIn, 2);
meanOut = NaN(1, nData);    stdOut = NaN(1, nData);

%% Compute the weighted mean and std
azmIn = u2u(azmIn, angleUnits, 'rad');
stdIn = u2u(stdIn, angleUnits, 'rad');

for iaz = 1:nData 
    % Find viable data entries
    notNaN = find(~isnan(azmIn(:,iaz))  &  ~isnan(stdIn(:,iaz)));     % find ~NaNs in the data

    if isempty(notNaN) == 0
        azmCmplx = exp(1i*azmIn(notNaN,iaz));           % Solve the problem in complex plane to  
                                                        % properly handle the 2*pi periodicity 
        [meanCmplx, stdCmplx] = wMean(azmCmplx, stdIn(notNaN,iaz));    % weighted std and mean
    else
        % stdIn = NaN for all data -- calculate the unweighted mean and std
        azmCmplx = exp(1i*azmIn(~isnan(azmIn(:,iaz)),iaz));        % ~NaN azimuths but NaN std's
        [meanCmplx, stdCmplx] = wMean(azmCmplx, NaN(size(azmCmplx)));  % unweighted std and mean
    end;    
       
    meanOut(1,iaz) = imag(log(meanCmplx));                                
    stdOut (1,iaz) = abs(stdCmplx);                                
    % This formulation ensures that, for example,
    % wAzim([2,2,2,4,5,5]', [1,1,1,1,1,1]', 'rad', []) = wAzim([2,4,5]', [1/3,1,1/2]', 'rad', []) and
    % wAzim([2,2,2,4,5,5]', [1,1,1,1,1,1]', 'deg', []) = wAzim([2,4,5]', [1/3,1,1/2]', 'deg', []),
    % however, clearly,
    % wAzim([2,4,5]', [1/3,1,1/2]', 'rad', []) ~= wAzim([2,4,5]', [1/3,1,1/2]', 'deg', [])
    
    if abs(meanCmplx) < tol
        fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
        fprintf(['>>> The data center of gravity is at %g from the origin', ...
                 ' of the unit circle \n'], abs(meanCmplx));
        fprintf('    Mean azimuth %g [%s] is ill-defined -- PAUSE(3) \n \n', ...
                u2u(meanOut(iaz,1), 'rad', angleUnits), angleUnits);   
        pause(3);
    end;        
end;

meanOut = u2u(meanOut, 'rad', angleUnits);
stdOut  = u2u( stdOut, 'rad', angleUnits);

end    % of the function