%% *xcorr1C*
% Cross-correlation or cross-sign of a single component data 

%%
% *Input:*

% traceData   - [header.noSample, header.noTrace] 3C trace array produced by function
%               'readSEGYSEG2'
% header      - structure containing the header information produced by function
%               'readSEGYSEG2'        
% noLag       - [scalar] the number of lags for cross-correlation 
%               (*) Matlab function 'xcorr' computes the positive and negative lags  
% iComp       - [scalar] data component to be cross-correlated     
% xcorrFlag   = [string] equal to 'Xcorr' or 'Xsign' to indicate whether the cross-correlation 
%               or cross-sign is to be computed  

%%
% *Output:*

% xcorData    - [header.noRec, header.noRec, 2*noLag+1] array of cross-correlations

%%
% *Author:* Vladimir Grechka 2012    

%%
function [xcorData] = xcorr1C(traceData, header, noLag, iComp, xcorrFlag)
%% Settings
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

xcorData = zeros(header.noRec, header.noRec, 1:2*noLag+1);

%% Extract 1C data from 3C shot gather
tic            
trace1C = traceData(:, iComp:3:header.noTrace);           
for itrace = 1:header.noRec
    tmp = trace1C(:,itrace);
    trace1C(:,itrace) = tmp/max(abs(tmp));      % normalize the traces before cross-correlation
end;
trace1C(isnan(trace1C)) = 0;                    % replace NaN trace values with zeros

%% Cross-correlate traces and write them into array 'xcorData'
for itrace = 1:header.noRec
    for jtrace = 1:header.noRec
        tmp1 = trace1C(:,jtrace); 
        tmp2 = trace1C(:,itrace);
        if strcmp(xcorrFlag, 'Xcorr') == 1 
            % do nothing
        elseif strcmp(xcorrFlag, 'Xsign') == 1 
            tmp1 = sign(tmp1); 
            tmp2 = sign(tmp2);
        else
            fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
            fprintf('>>> Parameter xcorrFlag = %s \n', xcorrFlag);
            fprintf('>>> Its supported values are ''Xcorr'' or ''Xsign'' \n \n');
              error('>>> STOP');
        end;
        
        xcorData(itrace, jtrace, 1:2*noLag+1) = xcorr(tmp1, tmp2, noLag, 'unbiased');
        
    end;
end;        
toc

end    % of the function

%%
