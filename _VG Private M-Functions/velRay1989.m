%% *velRay1989*
% Compute the wavefront normal and the phase and group velocities for a given ray direction with 
% the algorithm of Obolentseva and Grechka (1989) 

%%
% *Input:*

% Cij         - [6, 6] stiffness matrix
% r0          - [3, 1] unit ray vector specified by its spherical angles theta(i): 
%               r0 = [sin(theta(1))*cos(theta(2)), sin(theta(1))*sin(theta(2)), cos(theta(1))]' 
% n0          - [3, 1] initial guess for the wavefront normal n
%               If n0 is empty or n0 = NaN, r0 is taken as an initial guess for n
% waveType    - [scalar] equal to 1 (for the P-wave), or 2 (for the S1-wave), 
%               or 3 (for the S2-wave)
% tolParam    - [scalar] tolerance for finding r0
% noIter      - [scalar] number of iterations allowed to reach solution

%%
% *Output:*

% n           - [3, 1] wavefront normal specified by its spherical angles beta(i): 
%               n = [sin(beta(1))*cos(beta(2)), sin(beta(1))*sin(beta(2)), cos(beta(1))]' 
% V           - [scalar] phase velocity of wave waveType in direction n
% U           - [3, 1] polarization vector of the wave waveType
% g           - [3, 1] group-velocity vector of the wave waveType that has the ray direction r0
%               and the wavefront normal n  
% dg          - [2, 1] vector of derivatives d |g| / d theta(i)
% Kp          - [scalar] Gaussian curvature of the slowness surface
% misfit      - [scalar] norm(r - r0), where r is the obtained unit ray vector
% flagOut     - [scalar] convergence flag 
%               . flagOut = 1 - an accurate solution is found
%               . flagOut = 0 - no solution is found; 
%                               all other output parameters take their current values

%%
% *Author:* Vladimir Grechka 1989 2012 - 2014
%
% * Fortran version is published in Obolentseva and Grechka (1989)

%%
function [n, V, U, g, dg, Kp, misfit, flagOut] = velRay1989(Cij, r0, n0, waveType, tolParam, noIter)
%% Settings 
[~, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

tol1 = 1.e-6;   tol2 = 1.e-4;
iter = 0;       iterExtra = 0;    misfit = 1;     
shrinkOrExpand = [ [0.50; 2.0],   [1.0; 2.0],  [2.0; 2.0],  [1.0; 1.5], ...
                   [0.25; 1.0],   [0.5; 1.0],  [1.5; 1.0],  [2.0; 1.0], ...
                   [0.10; 0.5],   [0.5; 0.5],  [1.0; 0.5],  [2.0; 0.5], ...
                   [0.25; 0.25],  [1.0; 0.25], [0.1; 0.1],  [0.5; 0.1] ]; 

% Initial guess for the wavefront normal
if isempty(n0) == 1  ||  sum(isnan(n0)) > 0 
    n = r0;
else
    n = n0;
end;

% Remove singularity of the spherical coordinate system
if abs(r0(3)) > 1 - tol1;    r0(1) = tol2;   r0 = r0/norm(r0);    end;

%% Main iteration loop
while iter < noIter  &&  misfit > tolParam 
    iter = iter + 1;
    misfitOld = misfit;                             % save the value of variable misfit 

    % Compute the phase velocity
    [V1, U1] = velPhaseU(Cij, n);
    V = V1(waveType);   U = U1(:,waveType); 
    
    % Derivatives with respect to the pahse angles 
    % d n / d beta(i), d V / d beta(i), and d g / d beta(i) 
    [dn, spn] = dVecdAngle(n);
    [dVdbeta, dgdbeta] = dVGdBeta(Cij, n, V1, U1, waveType);        

    % Compute the group velocity (Grechka and Obolentseva, 1989, eq 7)
    g = V*n + dVdbeta(1)*dn(:,1) + dVdbeta(2)*dn(:,2)/spn^2;
    r = g/norm(g);                                  % new ray 
    if abs(r(3)) > 1 - tol1;    r(1) = tol2;   r = r/norm(r);   end;
    
    misfit = norm(r - r0);                          % update 'misfit'

    % Derivatives of |g| and g with respect to the ray angles
    % dg(i) = d |g| / d theta(i),  dgdtheta = d g / d theta(i)
    [dr, spr] = dVecdAngle(r);
    dg = -norm(g)*(dr'*n)/dot(r, n);                % eq 1.29 in Obolentseva and Grechka (1989)
    dgdtheta = r*dg' + norm(g)*dr;  
    
    % Derivatives dbdt(i,j) = d beta(i) / d theta(j)
    dbdt = pinv(dgdbeta)*dgdtheta;                  % eq 32 in Grechka and Obolentseva (1989)
    %eig(dot(r, n)*(dbdt)*spn/(dot(g, g)*spr))
    %det(dbdt)

    % Gaussian curvature of the group-velocity surface
    Kg = dot(r, n)*det(dbdt)*spn/(dot(g, g)*spr);   % eq 36 in Grechka and Obolentseva (1989)
    % Gaussian curvature of the slowness surface
    Kp = dot(r, n)^4/Kg;                            % eq 17 in Grechka and Obolentseva (1993)

    if Kp < 0                                       % negative Kp implies concavity of the 
        flagOut = 0;                                % slowness surface and, hence,   
        return;                                     % multivaluedness of the group-velocity 
    end;                                            % surface

    % Angular deviations dtheta1 and dtheta2 of the current ray r from the input ray r0
    dtheta1 = acos(r0(3)) - acos(r(3));
    var = (r(1)*r0(2) - r(2)*r0(1))/sqrt((1 - r(3)^2)*(1 - r0(3)^2));
    if abs(var) > 1;   var = (1 - tol1)*var/abs(var);    end;
    dtheta2 = asin(var);

    % Check the convergence
    if iter == 1  ||  misfit < misfitOld
        iterExtra = 0;
        dtheta1Old = dtheta1;                       % save angular deviations dtheta1 and
        dtheta2Old = dtheta2;                       % dtheta2
         nOld = n;                                  % save the wavefront normal and
        dnOld = dn;                                 % its derivatives

    else
        % Modify the step size over theta(i) in an attempt to get convergence
        iterExtra = iterExtra + 1;
        if iterExtra > size(shrinkOrExpand, 2)  ||  iter > noIter
            break;
        else
            dtheta1 = dtheta1Old*shrinkOrExpand(1,iterExtra);
            dtheta2 = dtheta2Old*shrinkOrExpand(2,iterExtra);
             n = nOld;
            dn = dnOld;
        end;
    end;
    
    % Update the wavefront normal
    nNew = n + (dn(:,1)*dbdt(1,1) + dn(:,2)*dbdt(2,1))*dtheta1 + ...
               (dn(:,1)*dbdt(1,2) + dn(:,2)*dbdt(2,2))*dtheta2;                
    nNew = nNew/norm(nNew);
    if abs(nNew(3)) > 1 - tol1;   nNew(1) = tol2;   nNew = nNew/norm(nNew);   end;
    n = nNew;
    
end; % of the while-loop

%% Set the output flag
if misfit < tolParam
    flagOut = 1;
else
    flagOut = 0; 
end

end    % of the function
