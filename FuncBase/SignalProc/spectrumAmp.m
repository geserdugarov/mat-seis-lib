%% *spectrumAmp*
% Calculating amplitude spectrum

%% 
% *Input:* 
% signal - [times (s), amplitudes] cutted windowed signal

%%
% *Output:*
% spAmp - [frequencies (Hz), amplitudes] amplitude spectrum

%%
% *Author:* Geser Dugarov 2016

%%
function [spAmp] = spectrumAmp(signal)

amplitudes = abs(fft(signal(:, 2)));
df = 1/(signal(end,1) - signal(1,1)); % frequency step
N  = round((size(signal,1)-1)/2); % number of frequencies from 0 to Nyquist frequency
amplitudes = amplitudes(1:N);
frequencies = df*(0:N-1)';
spAmp = [frequencies amplitudes];

end % of the function
