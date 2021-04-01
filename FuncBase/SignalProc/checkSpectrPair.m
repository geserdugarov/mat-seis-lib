%% *checkSpectrPair*
% Check length and discretization of two spectrums

%% 
% *Input:* 
% trace1 - [frequencies, amplitudes] spectrum
% trace2 - [frequencies, amplitudes] spectrum

%%
% *Output:*
% trace1 - [frequencies, amplitudes] spectrum
% trace2 - [frequencies, amplitudes] spectrum

%%
% *Author:* Geser Dugarov 2016

%%
function [trace1, trace2] = checkSpectrPair(trace1, trace2)

[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));

checkSignalAndDiscr(trace1);
checkSignalAndDiscr(trace2);

fmax1 = trace1(end,1);
fmax2 = trace2(end,1);
df1 = trace1(2,1)-trace1(1,1);
df2 = trace2(2,1)-trace2(1,1);
if (abs(fmax1-fmax2)/min([df1 df2])*100 > 0.1)
    if (100*abs(df1-df2)/mean([df1 df2]) < 0.1)
        fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
        fprintf('>>> WARNING: Not equal time length of signals. \n');
        fprintf('>>> WARNING: Skipping some points. \n \n');
        if (fmax1 > fmax2)
            temp = trace1(1:size(trace2,1),:);
            trace1 = temp;
        else
            temp = trace2(1:size(trace1,1),:);
            trace2 = temp;
        end
    else
        fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
        error('>>> Different discretization. Can not skip some points.');
    end
end

end % of the function
