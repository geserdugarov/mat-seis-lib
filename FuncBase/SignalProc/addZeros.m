%% *addZeros*
% Adding zeros in the end of the signal

%% 
% *Input:* 
% trace - [times (s), amplitudes] seismogram trace
% t     - [scalar (s)] time limit for adding zeros

%%
% *Output:*
% trace0 - [times (s), amplitudes] cutted windowed signal

%%
% *Author:* Geser Dugarov 2016

%%
function [trace0] = addZeros(trace, t)

trace0 = trace;
dt = trace(2,1)-trace(1,1);
if (trace(end,1)-trace(1,1) < t)
	num = round((trace(1,1)+t-trace(end,1)+1.0e-16)/dt);
    temp = trace(end,1);
    for i=1:num
        trace0 = [trace0; [temp+i*dt 0]];
    end
end

end % of the function
