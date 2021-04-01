%% *timePickManual*
% Manual time picking

%% 
% *Input:* 

% noRec        - [scalar] number of receivers in a seismogram
% missingValue - [scalar] indicating a missing value in the output array

%%
% *Output:*

% pickOut      - [noRec, 2] array whose columns contain the receiver numbers and the time picks 
%                made with Matlab function 'ginput' 

%%
% *Author:* Vladimir Grechka 2012

%% 
% *Assumptions and implications:*
%
% * the time axis of a seismogram is horizontal
% * the vertical axis is labeled in sequential receiver numbers rather than the receiver coordinates
% * time picks have the same time units as the traces from which the picks are made
% * *relaxing the assumptions above requires adding more input parameters, which I would like to 
%    avoid at this point*

%% 
% *Intended usage:*
%
% * display a seismogram 
% * select component and zoom into the area to be picked
% * run |timePickManual(noRec)| or |timePickManual(noRec, missingValue)|
% * make picks, hit ENTER when done
% * copy the picks into a spreadsheet for future use or save the output in a |.xlsx| or a |.mat| 
%   file
% * *the script below can be used to quickly estimate the apparent velocity*:
%   
%      t = timePickManual(noRec); t(any(isnan(t),2),:) = [], p = polyfit(dz*t(:,1), t(:,2), 1);
%      Vap = 1/p(1)

%%
function [pickOut] = timePickManual(noRec, missingValue)
%% Settings and defaults 
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
% narginTrue = nargin(thisFileName);   
% narginchk(narginTrue, narginTrue);

if nargin < 1  ||  isempty(noRec) == 1
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Missing or empty input parameter ''noRec'' \n \n');    
      error('>>> STOP');
end;

if nargin < 2  ||  isempty(missingValue) == 1
    missingValue = NaN;
end;
pickOut(:,1) = (1:noRec)';   
pickOut(:,2) = missingValue*ones(noRec, 1);

%% Time picking using 'ginput'
fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
fprintf('>>> Please make time picks \n');
tmp1 = ginput;
tmp2 = round(tmp1(:, 2));

%% Check the picks
mpicks = find(tmp2 > noRec);
if isempty(mpicks) == 0
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> WARNING: Time picks made at receivers \n');
    display('    rec #    time');
    for i = 1:length(mpicks)
        fprintf('     %g     %g \n', [tmp2(mpicks(i)), tmp1(mpicks(i), 1)']);
    end;
    fprintf('    exceeding noRec (= %g) will be discarded \n', noRec);    
    fprintf('>>> Please check the input -- PAUSE \n');    pause;
end;

%% Fill the output array
for irec = 1:noRec
    ipick = find(tmp2 == irec);
    if isempty(ipick) == 0;
        pickOut(irec,2) = tmp1(ipick(1), 1);
    end;
end;   

display('Picked times:');
display(pickOut(:,2));

end    % of the function
