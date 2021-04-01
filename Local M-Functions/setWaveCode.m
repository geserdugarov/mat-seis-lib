%% *setWaveCode*
% Specify wave codes 

%% 
% *NB:*
%
% * This script needs to be tailored to each specific data set

%%
% *Input:* none

%%
% *Output:*

% waveCode - [4, noWave] wave-code array
% noWave   - the number of wave codes
%            (*) Each column of 'waveCode' has four elements:
%               . waveCode(1,i) = 1, 2, or 3 - for P-, S1-, or S2-wave leaving a source 
%               . waveCode(2,i) = 0, 1, 2, or 3 - for reflected P-, S1-, or S2-wave
%               . waveCode(2,i) = 0 means no reflection so traveltime of a direct wave is computed                            
%               . waveCode(3,i) = 1, ..., noInt - the number of reflecting interface 
%               . waveCode(4,i) = 1, ..., noFault - the number of reflecting fault 
%            (*) If waveCode(2,i) = 0, the values of waveCode([3,4],i) are irrelevant    

%%
% *Author:* Vladimir Grechka 2012

%%
function [waveCode, noWave] = setWaveCode
%% Settings  
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

%% Input the 'waveCode'
global modelflag;
if modelflag==1
    waveCode = [1, 1, 2, 0;
                1, 2, 2, 0;
                1, 3, 2, 0]';
%     waveCode = [1, 2, 2, 0;
%                 1, 3, 2, 0]';
%     waveCode = [1, 2, 2, 0]';
elseif modelflag==2
    waveCode = [1, 1, 2, 0;
                1, 2, 2, 0;
                1, 3, 2, 0]';
elseif modelflag==3
    waveCode = [1, 2, 2, 0]';
%     waveCode = [1, 1, 2, 0;
%                 1, 2, 2, 0]';
elseif modelflag==4
%     waveCode = [1, 1, 1, 0;
%                 1, 2, 1, 0;
%                 1, 3, 1, 0]';
    waveCode = [1, 1, 1, 0]';
%     waveCode = [1, 2, 1, 0]';
elseif modelflag==5
    waveCode = [1, 1, 7, 0]';
elseif modelflag==6
    waveCode = [1, 1, 7, 0]';
elseif modelflag==7
    waveCode = [1, 1, 4, 0]';
end

% waveCode = [1, 0, 0, 0; ...             
%             2, 0, 0, 0; ...
%             1, 1, 4, 0; ...
%             2, 2, 4, 0]';  
noWave = size(waveCode, 2);

if size(waveCode, 1) ~= 4
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> The number of columns of array ''waveCode'' is incorrect \n');
    display(waveCode);
      error('>>> Correct the input');
end

end  % of the function