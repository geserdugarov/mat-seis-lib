%% *velPhaseU*
% Solve the Christoffel equation for the phase velocities and the polarization vectors 

%%
% *Input:*

% Cij - [6, 6] stiffness matrix
% n   - [3, 1] unit wavefront normal 
% tol - [scalar] tolerance for various numerical checks

%%
% *Output:*

% V   - [3, 1] phase velocities of isonormal plane waves 
%       (*) The velocities are sorted as V(P), V(SV), V(SH) for TI media and 
%           V(1) > V(2) >= V(3) for symmetries lower than TI 
% U   - [3, 3] polarization vectors (columns) of isonormal plane waves sorted the same way
%       as velocities 

%%
% *Author:* Vladimir Grechka 1998 2012 - 2014
%
% * Fortran version is published in Obolentseva and Grechka (1989)

%%
function [V, U] = velPhaseU(Cij, n, tol)
%% Settings and checks
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
% narginTrue = nargin(thisFileName);   
% narginchk(narginTrue, narginTrue);

if nargin < 3  ||  isempty(tol) == 1   
    tol = 1.e-12;                       % default for tol
end                   

% Check whether the wavefront normal is a unit vector
if abs(norm(n) - 1) > tol
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Wavefront normal n is not a unit vector within specified tolerance \n');
    fprintf('>>> abs(norm(n) - 1) = %g \n', abs(norm(n) - 1));
      error('>>> STOP');   
end

%% Phase-velocity calculation
if isISO(Cij, tol) == 1
    % Isotropy
    V(1) = sqrt(Cij(3,3));   V(2) = sqrt(Cij(5,5));   V(3) = V(2);
    U(:,1) = n;                                         % P-wave polarization  
    if abs(n(3)) > 1 - tol
        U(:,3) = [0; 1; 0];                             % SH-wave polarization  
        U(:,2) = [sign(n(3)); 0; 0];                    % SV-wave polarization    
    else
        U(:,3) = [-n(2); n(1); 0]/sqrt(1 - n(3)^2);     % SH-wave polarization
        U(:,2) = cross(U(:,1), U(:,3));                 % SV-wave polarization
    end

else
    % Symmetries lower than isotropy
    G = christoffel(Cij, n, n);                         % the Christoffel matrix
    [vec, V2] = eig(G);                                 % solve eigenvalue-eigenvector problem
    [V2sort, ind] = sort(diag(V2), 'descend');          % sort the roots so that
    V(1,:) = sqrt(V2sort);                              % V2(1) > V2(2) >= V2(3)
    U = vec(:, ind);
    
    if dot(U(:,1), n) < 0
        U(:,1) = -U(:,1);                               % flip the P-wave polarization if necessary 
        U(:,2) = -U(:,2);
    end

    [d1, spn, ~] = dVecdAngle(n);
    if dot(U(:,2), d1(:,1)) < 0
        U(:,2) = -U(:,2);                               % flip the S1-wave polarization if necessary 
        U(:,3) = -U(:,3);
    end

    if dot(U(:,3), d1(:,2)/spn) < 0
        U(:,3) = -U(:,3);                               % flip the S2-wave polarization if necessary 
        U(:,2) = -U(:,2);
    end

    %% Commented out for the Bakken
    % Find out whether a medium is TI and, if yes, possibly change the sorting of the shear-waves  
    [flag, tti2vti] = isTI(Cij, tol);
    if flag == 1    
        axisDir = tti2vti(:,3);                         % symmetry-axis vector in TI medium
        if abs(dot(U(:,2), axisDir)) < abs(dot(U(:,3), axisDir))
            % Switch the velocities and polarization vectors to ensure that the shear-waves 
            % are sorted as V(2) = V(SV), V(3) = V(SH) and U(:,2) = U(:,SV), U(:,3) = U(:,SH)
            tmpU = U(:,2);    U(:,2) = U(:,3);    U(:,3) = tmpU;
            tmpV = V(2);        V(2) = V(3);        V(3) = tmpV;
        end
    end
    %% End of the commented out block

end

%[n, U],  display('>>> Pause'); pause;

end    % of the function

%%
