%% *velRayVTI*
% Compute the wavefront normal and the phase and group velocities for a given ray direction in TI
% media

%%
% *Input:*

% Cij         - [6, 6] TTI stiffness matrix
% R           - [3, 3] rotation matrix that makes Cij VTI
% r           - [3, 1] ray vector 
% waveType    - [scalar] equal to 1 for the P-wave, or to 2 for the SV-wave, 
%               or to 3 (for the SH-wave)
% verboseFlag - [scalar] equal to 0 or 1 to controls the amount of information printed  
%               when the SV triplication is encountered

%%
% *Output:*

% n           - [3, 1] wavefront normal  
% V           - [scalar] phase velocity of wave waveType in direction n
% U           - [3, 1] polarization vector of the wave waveType
% g           - [3, 1] group-velocity vector of the wave waveType that has the ray direction r0
%               and the wavefront normal n  
% flagOut     - [scalar] output flag 
%               . flagOut = 1 - unique solution is found
%               . flagOut = 0 - input ray crosses triplication on the SV wavefront
%

%%
% *Author:* Vladimir Grechka 1998 2012 - 2014

%%
function [n, V, U, g, flagOut] = velRayVTI(Cij, R, r, waveType, verboseFlag)
%% Settings 
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

flagOut = 1;
r = r/norm(r);          % normalize r

%% Build rotation matrix R that makes Cij -- VTI and places r in the [x1, x3] plane
Rot = vec2axis(R(:,3), 3, r, 2);

%% Rotate TTI to VTI
C = bond(Cij, Rot');
r = Rot*r;

%% Grechka's (2012) solution for the wavefront normal, phase, and group velocities
if waveType == 3
    % SH-wave
    Psi = r(1)/r(3);                            % tangent of the ray angle
    Theta = Psi*C(5,5)/C(6,6);                  % tangent of the wavefront normal angle    
    n = [sin(atan(Theta)), 0, cos(atan(Theta))]';   % wavefront normal
    cosrn = dot(r, n);
    if cosrn < 0;    n = -n;    end;            % flip the sign of 'n' if necessary
    V = sqrt(C(6,6)*n(1)^2 + C(5,5)*n(3)^2);    % phase velocity
    g = r*V/abs(cosrn);                         % group-velocity vector
    U = [0, 1, 0]';                             % polarization vector

else
    [F, M] = setMatVTI(C(1,1), C(1,3), C(3,3), C(5,5));
    ani = cij2grechka(C); 
    [n, V, g, flagOut] = velVTI2012(F, M, [ani(1:3), ani(5)], r, waveType);
    if flagOut == 0
        U = NaN(3,1);
        if verboseFlag == 1
            fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
            fprintf('>>> Ray crosses triplication on the SV wavefront \n \n'); 
        end;              
%        error('>>> STOP');
    else
        [V1, U1] = velPhaseU(C, n);
        [~,pos] = min(abs(V1-V));
        U = U1(:,pos);                     % polarization vector
    end;
end;

%% Construct the output arrays
if flagOut == 1
    n = Rot'*n;    U = Rot'*U;    g = Rot'*g;
end;

end    % of the function