%% *spectrumPhase*
% Calculating phase spectrum

%% 
% *Input:* 
% signal - [times (s), amplitudes] cutted windowed signal

%%
% *Output:*
% spPhase - [frequencies (Hz), angles (deg)] phase spectrum

%%
% *Author:* Geser Dugarov 2016

%%
function [spPhase] = spectrumPhase(signal)

resFFT = fft(signal(:, 2));
angles = rad2deg(atan(-imag(resFFT)./real(resFFT)));
dt = signal(2,1) - signal(1,1); % time step
df = 1/(length(angles)*dt); % frequency step
N  = round(length(angles)/2); % number of frequencies from 0 to Nyquist frequency
angles = angles(1:N);
frequencies = df*(0:N-1)';
spPhase = [frequencies angles];

end % of the function
