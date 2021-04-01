%% *ratioSTALTA*
% Compute the ratio of amplitudes in short-to-long windows

%% 
% *Input:* 

% data         - [noSample, noTrace] data matrix containing 'noTrace' traces 'noSample' samples long
% traceArray   - [jTrace] array of traces on which computation is to be made; jTrace <= noTrace
% sampleCenter - [jTrace] array of centers (in samples) of the STA windows 
% windowSTA    - [scalar] half-window (in samples) around 'sampleCenter' for computing the 
%                average amplitude of the signal
% windowLTA    - [scalar] window (in samples) preceeding 'windowSTA' for computing the average 
%                amplitude of (presumably) noise

%%
% *Output:*

% SNR          - [1, noTrace] array of computed amplitude ratios

%%
% *Author:* Vladimir Grechka 2012 - 2014

%%
function SNR = ratioSTALTA(data, traceArray, sampleCenter, windowSTA, windowLTA)
%% Settings 
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

[noSample, noTrace] = size(data);
SNR = NaN(1, noTrace);

%% SNR computation
for itrace = traceArray
    s1 = sampleCenter(itrace) - windowSTA;
    s2 = sampleCenter(itrace) + windowSTA;
    s0 = s1 - windowLTA;
    if s0 < 1 || noSample < s2
        fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
        fprintf('>>> Violation of inequality 1 <= [sF - STA - LTA, sF - STA, sF, sF + STA] <= noSample: \n');
        fprintf('>>> 1 ?? [ %g, %g, %g, %g] ?? %g \n', s0, s1, sampleCenter(itrace), s2, noSample);
        beep on;  beep;  beep off;
          error('>>> STOP');
    end;

    tmpSTA = data(s1:s2, itrace);  
    tmpLTA = data(s0:s1, itrace);                                  
    SNR(itrace) = sqrt( dot(tmpSTA, tmpSTA)/length(tmpSTA) )/ ...
                  sqrt( dot(tmpLTA, tmpLTA)/length(tmpLTA) );
end;

%% Replace the remaining NaNs with ones
% SNR(isnan(SNR)) = 1;

end   % of the function 