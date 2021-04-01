%% *condNumHansen*
% Scale matrix in accordance with Hansen et al. (2012) to decrease its condition number

%%
% *Input:*

% Mat        - rectangular matrix whose condition number is to be minimized 
% printFlag  - [scalar] flag, equal to 0 or 1, to indicate verbose print 

%%
% *Output:*

% MatOut     - [size(Mat)] scaled version of input matrix 'Mat'
% scaleOut   - [1, size(Mat,2)] array of the inverse scaling factors for the parameter vector
% scaleFlag  - [scalar] flag equal to 0 when 'Mat' is numerically singular and 0 otherwise
% condNumber - [scalar] condition number of the scaled matrix 'Mat'

%%
% *Reference:*

% Hansen, P. C., V. Pereyra, and G. Scherer, 2012, Least Squares Data Fitting with Applications: 
% John Hopkins University Press.

%%
% *Note:*

% Algorithm of Hansen et al. is not guaranteed to always reduce the condition number. 
% Take, for example, matrix 
% a = [0.6463    0.2760    0.1626; ...
%      0.7094    0.6797    0.1190; ...
%      0.7547    0.6551    0.4984]
% Its cond(a) = 7.6799, whereas the output of 'condNumHansen' is condNumber = 8.5603.

%%
% *Author:* Vladimir Grechka 2014

%%
function [MatOut, scaleOut, scaleFlag, condNumber] = condNumHansen(Mat, printFlag)
%% Settings and defaults
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

tolMax    = 1.e+8/norm(Mat);  % tolerance intended to remove the amplification of round-off errors
scaleFlag = 1;
scaleOut  = ones(1,size(Mat,2));

if printFlag == 1
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    condNumber = cond(Mat);
    fprintf('>>> Condition number prior to scaling = %g \n', condNumber);
end;

%% Apply scaling proposed in Theorem 48 in Hansen et al. (2012, p. 56) 
for i = 1:size(Mat,2)
    scaleOut(1,i) = 1/norm(Mat(:,i));
    if abs(scaleOut(1,i)) > tolMax
        scaleFlag = 0;
        scaleOut(1,i) = tolMax;
        fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
        display('>>> WARNING:');
        fprintf('        Data do not constrain variable %g \n', i);
        display('>>> RECOMMENDED ACTION:');
        display('        (1) Check whether the issue is related to different scaling of columns');
        display('        (2) If yes, ignore the message; if no, stop the inversion and');
        fprintf('            rerun it after removing parameter %g from the unknowns \n', i);
%        fprintf('>>> PAUSE \n \n');  pause;
    end;
end;
MatOut = Mat*diag(scaleOut);  condNumber = cond(MatOut);  

if printFlag == 1
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Condition number after scaling = %g \n', condNumber);
end;

end    % of the function
