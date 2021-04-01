%% *calcTracesZ*
% Prepare data arrays for 3C seismograms from ray data

%%
% *Input:*

% t    - [1, noRay]       array of arrival times, s 
% U    - [3, noRay]       polarization vectors
% R    - [noLayer, noRay] reflection/transmission coefficients
% L    - [1, noRay]       geometric spreading coefficients
% type - [string]         type of wavelet
% T    - [scalar]         seismogram length, s
% dT   - [scalar]         seismogram discretization, s

%%
% *Output:*

% x   - [1,T/dT] time array
% y   - [1,T/dT] array of wavelet values

%%
% *Author:* Geser Dugarov 2019

function [time, Z] = calcTracesZ(t, U, R, L, type, T, dT)

noRay = numel(t);
noLayer = size(R,1);
time = (0:dT:T)';
Z = zeros(numel(time),noRay);

for ray = 1:noRay
    % constructing wavelet
    if strcmp(type,'ricker')
        Twind = 0.1; % length, s
        f = 25; % frequency, Hz
        tau = rem(t(ray),dT);
        x = -Twind/2:dT:Twind/2;
        y = (1.0 - 2.0*pi^2*f^2.*(x-tau).^2).*exp(-pi^2*f^2.*(x-tau).^2);
    end
    
    % accounting for reflection/transmission
    for layer = 1:noLayer
        y = y*real(R(layer,ray)); % didn't count phase shifting !!!
    end
    % accounting for geometric spreading (not realized)
    y = y/L(1,ray);
    
    % adding wavelet to seismogram
    pos = floor(t(ray)/dT)+1; % +1 for 0 accounting
    Z(pos+1:pos+numel(x),ray) =  y*U(3,ray);
end

end

%