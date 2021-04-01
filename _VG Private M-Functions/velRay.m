%% *velRay*
% Compute the wavefront normal and the phase and group velocities for a given ray direction 

%%
% *Input:*

% Cij         - [6, 6] stiffness matrix
% r0          - [3, 1] unit ray vector specified by its spherical angles theta(i): 
%               r0 = [sin(theta(1))*cos(theta(2)), sin(theta(1))*sin(theta(2)), cos(theta(1))]' 
% waveType    - [scalar] equal to 1 (for the P-wave), or to 2 (for the S1-wave), 
%               or to 3 (for the S2-wave)
% tolParam    - [scalar] tolerance for finding r0
% noIter      - [scalar] number of iterations allowed to reach solution
% verboseFlag - [scalar] equal to 1 or 0 that controls the amount of printed information 
%               describing a poor convergence

%%
% *Output:*

% n           - [3, 1] wavefront normal  
% V           - [scalar] phase velocity of wave waveType in direction n
% U           - [3, 1] polarization vector of the wave waveType
% g           - [3, 1] group-velocity vector of the wave waveType that has the ray direction r0
%               and the wavefront normal n  
% dg          - [2, 1] vector of derivatives d |g| / d theta(i)
% Kp          - [scalar] Gaussian curvature of the slowness surface
% flagOut     - [scalar] convergence flag 
%               . flagOut = 1 - an accurate solution is found
%               . flagOut = 0 - no solution is found; 
%                               all other output parameters take their current values
%
% (*) An alternative to the algorithm implemented here might be the minimization of |r0 - r(n)| in
%     terms of the wavefront normal angles beta. The search can be facilitated by the derivatives
%     dg/dbeta computed in function 'dVGdBeta'

%%
% *Author:* Vladimir Grechka 1998 2012 - 2014

%%
% *Comment:* 
%
% * An alternative to the algorithm implemented here might be the minimization of |norm(r0 - r(n))|
%   in terms of the wavefront normal angles |beta|. The search can be facilitated by the gradient
%   |dg/dbeta| computed in function |dVGdBeta|

%%
function [n, V, U, g, dg, Kp, flagOut] = velRay(Cij, r0, waveType, tolParam, noIter, verboseFlag)
%% Settings 
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

flagOut = 1;   tol = 1.e-6;    Kp = NaN;
r0 = r0/norm(r0);                               % normalize r0

%% Calculate the wavefront normal and velocities
flagISO = isISO(Cij, tol);

if flagISO == 1                                 % isotropy
    n = r0;
    [V1, U1] = velPhaseU(Cij, n);              
    V = V1(waveType);   U = U1(:,waveType);    g = V*r0;    dg = [0, 0]';    Kp = V^2;
   
else
    %% Comment out this block and also "end" on line 129 to enforce the S1-S2 mode selection for TI
    [flagVTI, RotVTI] = isTI(Cij, tol);
    if flagVTI == 1                             % transverse isotropy
        [n, V, U, g, flagOut] = velRayVTI(Cij, RotVTI, r0, waveType, verboseFlag);
        r = r0;
        if flagOut == 1
            [dr, ~] = dVecdAngle(r);
            dg = -norm(g)*(dr'*n)/dot(r, n);
        else
            dg = NaN(2,1);
        end;
    else                                        % symmetries lower than TI
    %% End of the commented out block

        iter = 0;   misfit = 1;                 % set parameters for simple iterartion 
        Rot = vec2axis(r0, 1, r0, 2);           % make r0 = [1, 0, 0]' to avoid the degeneracy          
        r0 = Rot*r0;                            % of spherical coordinates at the vertical
        [CijRot, ~] = bond(Cij, Rot');  
        n = r0;                                 % initial guess for the wavefront normal

        % Iterartions to find the wavefront normal n
        while iter < noIter  &&  misfit > tolParam 
            iter = iter + 1;
            misfit0 = misfit;

            % Compute the phase velocity
            [V1, U1] = velPhaseU(CijRot, n);

            % Compute the group velocity
            g1 = velGroup(CijRot, n, V1, U1, waveType);
            V = V1(waveType);   U = U1(:,waveType);   g = g1(:,waveType);   r = g/norm(g);
            misfit = norm(r - r0);
            % Check the convergence
            if misfit > 0.8*misfit0 
                % Iterations converge slowly or diverge --> switch to the 1989 algorithm
                n0 = n;
                [n, V, U, g, dg, Kp, ~, flagOut] = ...
                    velRay1989(CijRot, r0, n0, waveType, tolParam, noIter);
                if flagOut == 0   
                    if verboseFlag == 1
                        fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, ...
                                                                    [thisFileName, '.m']));
                        fprintf('>>> Iterations have not converged \n'); 
                        if Kp < 0
                            fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, ...
                                                                        [thisFileName, '.m']));
                            fprintf('>>> Slowness surface is concave: \n');
                            fprintf('    Its Gaussian curvature = %g \n', Kp); 
                        end;
                    end;
                    return;
                end;

            else
                % Iterations converge rapidly --> proceed
                rotAxis = cross(r, r0);
                rotAngl = asin(norm(rotAxis));
                if norm(rotAxis) > tolParam
                    R = rotateAroundAxis(rotAxis/norm(rotAxis), rotAngl, 'rad');
                else 
                    R = eye(3);
                end;
                n = R*n;                % rotate n with rotation matrix R that turns r into r0
            end;
       end;

        if flagOut == 1
            % Rotate output quantities into the original coordinate frame
            n = Rot'*n;    U = Rot'*U;    g = Rot'*g;    r = g/norm(g);
        end;
    end;

    %% Calculate derivative of the group velocity with respect to the ray angles
    if flagOut == 1  
        dr = dVecdAngle(r);                     % dr(:,i) = d r(:)/d theta(i)
        dg = -norm(g)*(dr'*n)/dot(r, n);        % eq 1.29 in Obolentseva and Grechka (1989)
        if isnan(Kp) == 1;
            [~, ~, ~, ~, ~, Kp, ~, ~] = velRay1989(Cij, r, n, waveType, tolParam, 1);
        end;
    end;
end;

end    % of the function
