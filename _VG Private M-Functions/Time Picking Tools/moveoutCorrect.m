%% *moveoutCorrect*
% Moveout-correct seismic gather

%% 
% *Input:*

% seismicGather    - [:, noTrace] seismic data matrix containing 'noTrace' seismic traces 
% moveout          - [noTrace] array of moveout corrections (in samples)
% moveoutNew       - [noTrace] array containing the new moveout (in samples) onto 
%                    which 'moveout' is mapped

%% 
% *Output:*

% movCorrGather     - [:, noTrace] moveout-corrected 'seismicGather'
% moveoutCorrection - [noTrace] array of the applied moveout corrections (in samples)

%%
% *Author:* Vladimir Grechka 2013 2014

%%
function [movCorrGather, moveoutCorrection] = moveoutCorrect(seismicGather, moveout, moveoutNew)
%% Settings and defaults
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

noSample = size(seismicGather, 1);  noTrace = size(seismicGather, 2);  
movCorrGather = zeros(noSample, noTrace);
%movCorrGather = seismicGather;

if any(moveoutNew) < 0 || any(moveoutNew) > noSample
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Sample(s) of ''moveoutNew'' are out of range [1, %g] \n', noSample);
    display(moveoutNew);
      error('>>> STOP');
end;
    
%% Calculate the time shifts and apply the moveout correction
moveoutCorrection = moveoutNew - moveout;
for itrace = 1:noTrace
    if isnan(moveoutCorrection(itrace)) == 0
        if moveoutCorrection(itrace) >= 0
            movCorrGather(1 + moveoutCorrection(itrace) : noSample, itrace) = ...
            seismicGather(1 : noSample - moveoutCorrection(itrace), itrace);
        else
            movCorrGather(1 : noSample + moveoutCorrection(itrace), itrace) = ...
            seismicGather(1 - moveoutCorrection(itrace) : noSample, itrace);
        end;
    else
        movCorrGather(:, itrace) = seismicGather(:, itrace);  % no moveout correction for NaNs
    end;
end;        

end  % of the script
