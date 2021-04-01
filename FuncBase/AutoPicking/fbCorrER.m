%% *fbCorrER*
% Correcting first brakes using ER or dER functions

%% 
% *Input:* 
% ER - [times (s), ER values]
% fb - [scalar (s)] approximate first brake time
% wind - [scalar (s)] ER/dER time window
% type - [scalar] correction type, 1 - ER, 2 - dER

%%
% *Output:*
% fbnew - [scalar (s)] updated first brake time

%%
% *Author:* Geser Dugarov 2016

%%
function [fbnew] = fbCorrER(ER, fb, wind, type)
    valuenum = size(ER,1);
    [~,fbpos] = min(abs(ER(:,1) - fb)); % index for first brake
    winddiscr = round(wind/(ER(2,1)-ER(1,1))); % convert time to discretes
    switch type
        case 1
            % correcting using ER
            if(fbpos-winddiscr > 0)
                bounddata = ER(fbpos-winddiscr:fbpos-floor(winddiscr/2),2);
                bound = mean(bounddata)+3*std(bounddata);
                bounddata = ER(fbpos-winddiscr:fbpos,2);
                pos = find(abs(bounddata)>bound);
                if isempty(pos)
                    fbnew = ER(fbpos,1);
                else
                    fbnew = ER(fbpos-winddiscr+pos(1,1),1);
                end
                clear bound bounddata pos;
            else
                fbnew = ER(fbpos,1);
            end
        case 2
            % correcting using dER
            if(fbpos-winddiscr > 0)
                dER = circshift(ER(:,2),-1)-ER(:,2);
                bounddata = dER(max([1,fbpos-winddiscr]):min([fbpos+winddiscr,valuenum]),1);
                bound = max(bounddata)*0.1;
                pos = find(abs(bounddata)>bound);
                if isempty(pos)
                    fbnew = ER(fbpos,1);
                else
                    fbnew = ER(fbpos-winddiscr+pos(1,1),1);
                end
                clear bound bounddata pos;
            else
                fbnew = ER(fbpos,1);
            end
        otherwise
            error('Wrong ER correction type: %s \n',num2str(type));
    end
end % of the function