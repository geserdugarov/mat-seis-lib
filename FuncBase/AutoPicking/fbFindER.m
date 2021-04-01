%% *fbFindER*
% Correcting first brakes using ER or dER functions

%% 
% *Input:* 
% ER - [times (s), ER values]
% wind - [scalar (s)] time window for std and mean
% skipSmall - true - skip first small signals
% windSkip - [scalar (s)] time window for searching max amplitude for 
%                         skipping first small signals
% trace - [times (s), amplitudes] seismogram
% badpos - [times (s)] bad first brakes

%%
% *Output:*
% fb - [scalar (s)] first brake time

%%
% *Author:* Geser Dugarov 2016

%%
function [fb, EPS] = fbFindER(ER, wind, skipSmall, windSkip, trace, badpos)
    valuenum = size(ER,1);
    winddiscr = round(wind/(ER(2,1)-ER(1,1))); % convert time to discretes
    windSkipD = round(windSkip/(ER(2,1)-ER(1,1)));
    badpos(isnan(badpos)) = [];
    if ~isempty(badpos)
        badposD = round(badpos/(ER(2,1)-ER(1,1)));
    end
    stdNe = ceil(winddiscr/2)-1; % movingstd window length (~ n/2)
    tracedata = ER(:,2);
    stds = movingstd(tracedata,stdNe,'central'); % from center to left and right
    means = movingmean(tracedata,winddiscr); % from left to right (if even then changed to n-1)
    stdsMat  = NaN(valuenum,2*stdNe+1);
    meansMat = NaN(valuenum,2*stdNe+1);
    stdsMat(:,stdNe+1)  = stds;
    meansMat(:,stdNe+1) = means;
    % search min element position k in sliding std window and
    % taking k-th element in mean window
    for slide=1:stdNe
        stdsMat(:,stdNe+1-slide) = circshift(stds,-slide);
        stdsMat(end-stdNe+1:end,stdNe+1-slide) = NaN;
        stdsMat(:,stdNe+1+slide) = circshift(stds,slide);
        stdsMat(1:stdNe,stdNe+1+slide) = NaN;
        meansMat(:,stdNe+1-slide) = circshift(means,-slide);
        meansMat(end-stdNe+1:end,stdNe+1-slide) = NaN;
        meansMat(:,stdNe+1+slide) = circshift(means,slide);
        meansMat(1:stdNe,stdNe+1+slide) = NaN;
    end
    clear slide;
    [~, stdind] = min(stdsMat, [], 2);
    EPS = diag(meansMat(:,stdind));
    clear stdind;
    deps = circshift(EPS(:,1),-1)-EPS(:,1);
    deps(end,1) = NaN;
    
    if ~skipSmall
        [~,pos] = max(deps);
        if ~isempty(badpos)
            if min(abs(badposD - pos(1,1))) > winddiscr/2
                fb = ER(pos(1,1),1);
                return;
            end
        end
    end

    deps(isnan(deps)) = 0;
    checkposnum = max([3 numel(badpos)+5]); ampl = NaN(checkposnum,1);
    [val,pos] = sort(deps,'descend'); pos = pos(1:checkposnum,1); val = val(1:checkposnum,1);
    temp = EPS(max([ones(checkposnum,1) pos-max([1,round(windSkipD/4)])],[],2),1);
    epsDiffTo0 = (EPS(pos+1,1)-temp)./max([abs(temp),repmat(1.0e-6,checkposnum,1)],[],2);
    for i = 1:checkposnum
        ind = max([1,pos(i,1)-windSkipD]):min([pos(i,1)+windSkipD,valuenum]);
        ampl(i,1) = max(abs(trace(ind,2)));
        % skip all bad first breaks
        if ~isempty(badpos)
            if min(abs(badposD - pos(i,1))) < winddiscr
                ampl(i,1) = 0;
                epsDiffTo0(i,1) = 0;
            end
        end
    end
    % skip all with EPS not close to 0
    [~,ind] = max(epsDiffTo0);
    for i = 1:checkposnum
        if epsDiffTo0(i,1) < epsDiffTo0(ind,1)*0.1;
            ampl(i,1) = 0;
        end
    end
    % skip all with small dEPS
    ind = find(ampl(:,1),1,'first');
    for i = checkposnum:-1:1
        if val(i,1) < val(ind,1)*0.4;
            ampl(i,1) = 0;
        end
    end
    % skip all with small amplitudes
    [~,ind] = max(ampl);
    for i = 1:checkposnum
        if ampl(i,1) < ampl(ind,1)*0.4;
            ampl(i,1) = 0;
        end
    end

    pos2 = find(ampl(:,1),1,'first');
    fb = ER(pos(pos2,1),1);

end % of the function