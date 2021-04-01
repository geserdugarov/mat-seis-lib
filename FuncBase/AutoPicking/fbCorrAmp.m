%% *fbCorrAmp*
% Correcting first brakes using amplitudes

%% 
% *Input:* 
% trace - [times (s), amplitudes]
% fb - [scalar (s)] approximate first brake time
% wind - [scalar (s)] time window

%%
% *Output:*
% fbnew - [scalar (s)] updated first brake time

%%
% *Author:* Geser Dugarov 2016

%%
function [fbnew] = fbCorrAmp(trace, fb, wind)
    valuenum = size(trace,1);
    [~,fbpos] = min(abs(trace(:,1) - fb)); % index for first brake
    winddiscr = round(wind/(trace(2,1)-trace(1,1))/2); % convert time to discretes
    if(fbpos-winddiscr > 0 && fbpos+winddiscr < valuenum)
        cutSig = trace(fbpos-winddiscr:fbpos+winddiscr,:);
        intpoly = polyfit(cutSig(:,1),cutSig(:,2),1);
        if ~(intpoly(1,1) == 0)
            fbnew = -intpoly(1,2)/intpoly(1,1);
        else
            fbnew = fb;
        end
    else
        fbnew = fb;
    end

end % of the function