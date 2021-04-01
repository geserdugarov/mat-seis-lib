%% *QfromSpRatio*
% Estimation inverse quality factor using spectral ratio method

%% 
% *Input:* 
% spAmpBef - [frequencies (Hz), amplitudes] amplitude spectrum of reference
%                                           signal
% spAmpAft - [frequencies (Hz), amplitudes] amplitude spectrum of signal
%                                           after attenuation
% frRange  - [f1(Hz) f2(Hz)] frequency range for linear interpolation
% t        - [scalar (s)] signal travel time from source to receiver
% flag     - [on off] plot picture
% figID    - figure ID
% figpos   - figure position
% xrange   - plotting range for X axis

%%
% *Output:*
% invQ    - [scalar] inverse quality factor
% invQerr - [scalar] inverse quality factor estimation error
% LnSpOut - [frequencies (Hz), amplitudes] logarithm of the spectral ratio
% linfit  - [scalar scalar] coefficients of linear regression for LnSpOut

%%
% *Author:* Geser Dugarov 2016

%%
function [invQ, invQerr, LnSpOut, linfit] = QfromSpRatio(spAmpBef, spAmpAft, frRange, t, flag, figID, figpos, xrange)

[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));

if (frRange(1)>frRange(2) || frRange(1)<0)
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    error('>>> Bad frequency range. \n \n');
end

freq = spAmpBef(:,1);
LnSp1 = log(spAmpBef(:,2));
LnSp2 = log(spAmpAft(:,2));
freqLn = [];
LnSp1Sp2 = [];
for i=1:size(freq,1)
    if (frRange(1)<=freq(i,1) && freq(i,1)<=frRange(2))
        freqLn = [freqLn; freq(i,1)];
        LnSp1Sp2 = [LnSp1Sp2; LnSp1(i,1)-LnSp2(i,1)];
    end
end
ws = warning('off','all');
linfit = polyfit(freqLn,LnSp1Sp2,1);
warning(ws);
k = linfit(1);
invQ = k/pi/t;
LnSpOut = [freq LnSp1-LnSp2];

% error estimation
n = size(freqLn,1);
disp = norm(polyval(linfit,freqLn)-LnSp1Sp2)^2/(n-2);
kerr = n*disp/(n*norm(freqLn)^2-sum(freqLn)^2);
invQerr = kerr/pi/t;

if strcmp(flag,'on')
    if (isempty(figID) && isempty(figpos))
        figure;
        hold on;
        title(['Q^{-1} = ',num2str(invQ)], 'Interpreter', 'tex')
    else
        setFigure(figID, 2, figpos, {['Q^{-1} = ',num2str(invQ)]}, ...
                'Frequencies, Hz', '', [], 16, 'tex');
    end
    line1 = plot(freq,LnSp1,'b');
    line2 = plot(freq,LnSp2,'r');
    line3 = plot(freq,LnSp1-LnSp2,'-*g');
    plot(frRange,polyval(linfit,frRange),'k','LineWidth',1.5);
    if isempty(xrange)
        xlim([0 3*frRange(2)]);
    else
        xlim(xrange);
    end
    set(gca, 'Ydir', 'normal');
    hold off;
    grid on;
    legend([line1, line2, line3], {'ln(SpRef)','ln(SpAft)','ln(SpR/SpA)'}, ...
           'Location','southeast','Orientation','horizontal');
end

end % of the function