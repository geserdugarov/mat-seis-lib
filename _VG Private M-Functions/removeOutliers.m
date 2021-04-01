%% *removeOutliers*
% Remove outliers from a vector 

%%
% *Input:*

% vecIn    - input vector
% level    - [scalar] level measured in standard deviations of 'vecIn' above which an element
%            is defined as an outlier and replaced with NaN 

%%
% *Output:*

% vecOut   - [size(vecIn)] vector in which the outliers are replaced with NaNs 
% outliers - [size(vecIn)] vector containing NaNs at the positions of outliers in 'vecIn'; 
%            the rest of 'outliers' is filled with zeros

%%
% *Author:* Vladimir Grechka 2014

%%
function [vecOut, outliers] = removeOutliers(vecIn, level)
%% Settings
[~, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

vecOut = vecIn;
outliers = zeros(size(vecIn));

meanVec = mean(vecIn(isnan(vecIn) == 0));
stdVec  =  std(vecIn(isnan(vecIn) == 0));

%% Identify outliers and replace them with NaNs
for i = 1:length(vecOut)
    if isnan(vecOut(i)) == 1 || abs(vecOut(i) - meanVec) > level*stdVec
        outliers(i) = NaN;
        vecOut(i)   = NaN;
    end;
end;

end  % of the function