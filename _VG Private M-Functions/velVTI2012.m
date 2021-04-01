%% *velVTI2012*
% Calculate the wavefront normals, phase, and group velocities for a given ray in VTI medium 

%% 
% *Input:*

% F        - [3, 3] matrix computed in function 'setMatVTI' that represents the Christoffel equation 
%            in the form [1, p1^2, p1^4] * F * [1, p3^2, p3^4]' = 0, where p = [p1, 0, p3] is the 
%            slowness vector
% M        - [7, 3] matrix computed in function 'setMatVTI' that relates Psi = tan(group angle) to 
%            Theta = tan(phase angle) as 
%            [Theta^0, Theta^1, ... Theta^6] * M * [Psi^0, Psi^1, Psi^2]' = 0
% ani      - [1, 4] vector containing Thomsen parameters [Vp0, Vs0, epsilon, delta] 
% r        - [3, 1] ray direction 
% waveType - [scalar] equal to 1 for the P-wave or to 2 for the SV-wave

%% 
% *Output:*

% n        - [3, 1] wavefront normal 
% V        - [scalar] phase velocity
% g        - [3, 1] group-velocity vector
% flagOut  - [scalar] output flag 
%            . flagOut = 1 - unique solution is found
%            . flagOut = 0 - input ray crosses triplication on the SV wavefront

%% 
% *Comment:*
%
% * Function |velVTI| is a modification of function with the same name published in Grechka (2012)
% that handles the P- and SV-waves only

%%
% *Author:* Vladimir Grechka 2012 - 2014

%%
function [n, V, g, flagOut] = velVTI2012(F, M, ani, r, waveType)
%% Settings
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

if waveType ~= 1  &&  waveType ~= 2
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Variable waveType = %g. It''s proper values are 1 or 2 \n \n', waveType); 
      error('>>> STOP');
end;

Theta = NaN(1,4);   Vph = NaN(1,4);
tol = 5.e-4;                    % parameter for identifying real-valued and double roots in 
                                % [Theta^0, Theta^1, ... Theta^6] * M * [Psi^0, Psi^1, Psi^2]' = 0
% !! The value tol = 5.e-4 is the same as in Geophysics G2P paper. Smaller values of tol, e.g.,
% !! tol = 1.e-6 cause missing rela-valued roots 'realRoots' on line 63 
r   = r/norm(r);                % normalize r
r(abs(r) < eps)   = eps;        % perturb variable 'r' to avoid division by zero
r(abs(r) > 1/eps) = 1/eps;
Psi = r(1)/r(3);                % tangent of the ray angle

%% Compute Theta-coefficients of the polynomial [Theta^0, Theta^1, ... Theta^6] * M * [Psi^0, Psi^1, Psi^2]' = 0
pol = NaN(1,0);    powers = NaN(1,3);
for i=1:3;   powers(1,i) = Psi^(i-1);   end;
for i=1:7;   pol(i) = dot(M(i,:), powers, 2);   end;

%% Solve equation [Theta^0, Theta^1, ... Theta^6] * M * [Psi^0, Psi^1, Psi^2]' = 0 for Theta 
allRoots = roots(fliplr(pol));
realRoots = real(allRoots(abs(imag(allRoots)) < tol));      % real-valued roots

% Find double roots for which the Christoffel equation is to be solved
doubleRoots = NaN(size(realRoots));
for i=1:length(realRoots)-1
    for j=i+1:length(realRoots)
        if abs(realRoots(i) - realRoots(j)) < tol*(abs(realRoots(i)) + abs(realRoots(j)))
            doubleRoots(i) = i;
            doubleRoots(j) = j;
        end;
    end;
end;

%% Compute the P- and SV-wave phase velocities
noWave = 0;     doubleRootIndex = 0;
for i=1:length(realRoots)
    noWave = noWave + 1;
    Theta(noWave) = realRoots(i);

    if isnan(doubleRoots(i)) == 1
        % Single root: compute the phase velocity as V = V(Theta, Psi) 
        p32 = (F(2,1)*Theta(noWave) - F(1,2)*Psi)/ ...         % squared vertical slowness
              (2*F(1,3)*Psi - F(2,2)*Theta(noWave) + ...
               F(2,2)*Psi*Theta(noWave)^2 - 2*F(3,1)*Theta(noWave)^3);
        Vph(noWave) = 1/sqrt(p32*(1 + Theta(noWave)^2));       % phase velocity
    else
        % Double root: compute the phase velocity from the Christoffel equation
        if abs(ani(3)) < tol  &&  abs(ani(4)) < tol  
            Vph(noWave) = ani(noWave);            % isotropy
        else
            f = 1 - (ani(2)/ani(1))^2;            % VTI
            tmp = atan(Theta(noWave));    s2 = sin(tmp)^2;    c2 = cos(tmp)^2;
            Vph(noWave) = ani(1)*sqrt(1 + ani(3)*s2 - f/2 + ...
                (-1)^doubleRootIndex*(f/2)*sqrt(1 + (4*s2/f)*(2*ani(4)*c2 - ani(3)*(c2 - s2)) + ...
                                                4*ani(3)^2*s2^2/f^2));
            doubleRootIndex = doubleRootIndex + 1;
        end;        
    end;
end;            

Vph = Vph(~isnan(Vph));       % remove NaN's

%% Sort the roots to make the first one to correspond to the P-wave
[Vph, index] = sort(Vph, 'descend');

%% Compute the rays, the wavefront normals, and the group velocities
if waveType == 2  &&  length(Vph) > 2
    flagOut = 0;        % SV-ray crossing a triplication
    n = NaN(3,1);    V = NaN;    g = NaN(3,1);
    return;
else
    flagOut = 1;        % Unique solution is found
    tmp = atan(Theta(index(waveType)));
    n = [sin(tmp), 0, cos(tmp)]';
    cosrn = dot(r, n);
    if cosrn < 0;    
        n = -n;                                 % flip the sign of 'n' if necessary
    end;            
    V = Vph(waveType);
    g = r*V/abs(cosrn);                         % group-velocity vector
end;

end  % of the function