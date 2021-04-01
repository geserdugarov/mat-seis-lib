%% *nmoCylInt*
% Compute an interval NMO cylinder (Grechka and Tsvankin, 2002) 

%%
% *Input:*
% Cij      - [6, 6] stiffness matrix in Voigt notation
% p        - [3, 1] slowness vector
% waveType - [scalar] equal to 1 (for the P-wave), 2 (for the S1-wave or SV-wave), 
%            or 3 (for the S2-wave or SH-wave) 

%%
% *Output:*
% Ucyl     - [3, 3] symmetric matrix representing an interval NMO cylinder
% q        - [2, 1] array of derivatives of the vertical slowness dp(3)/dp(i),  (i = 1, 2)

%%
% *Author:* Vladimir Grechka 1998, 2014

%%
% *Comments:*
%
% * Original Matlab version, called 'intCyl', is a part of the 'ART' package freely  
%   distributed by the CWP, CSM
%
% * Input parameter 'waveType' is only used to discriminate the shear modes at a singularity

%%
function [Ucyl, q] = nmoCylInt(Cij, p, waveType)
%% Settings  
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

% global G dG d2G

tol = 1.e-12;  tolSWS = 1.e-4;
q = zeros(2,1);  f = zeros(3,1);  ff = zeros(3,3);  Q = zeros(2,2);  

if isISO(Cij, tol) == 1
    %% Isotropy
    Ucyl = dot(p, p)*eye(3) - kron(p, p');
    for i = 1:2;
        q(i) = -p(i)/p(3);   
    end;
    
else
    %% Anisotropy
    % The Christoffel matrix in slowness G(i,j) = Cij([ik],[mj])*p(k)*p(m),
    G = christoffel(Cij, p, p) - eye(3);  

    % Find out if the slowness vector 'p' corresponds to a shear-wave singularity
    eigG = sort(abs(eig(G)), 'ascend');     % abs(eigG(1)) is always zero because the slowness
    if eigG(2) < tol                        % vector p satisfies the Christoffel equation
        fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
        fprintf(['>>> WARNING: Shear-wave singularity is encountered --> Small directional \n', ...
                 '             perturbation is applied to move away from the singularity \n \n']);
        wfn = p/norm(p);
        wfn = wfn + tolSWS*randn(size(p));  % perturb the wavefront normal
        wfn = wfn/norm(wfn);
        [Vphase, ~] = velPhaseU(Cij, wfn, tol);
        p = wfn/Vphase(waveType);           % get the new slowness vector
        G = christoffel(Cij, p, p) - eye(3);  
    end;
    
    % Derivatives of the Christoffel matrix G
    % dG(i,j,k) = [d G(i,j)]/[d p(k)], (i, j, k = 1, ..., 3)
    % d2G(i,j,k,m) = [d^2 G(i,j)]/[d p(k) d p(m)], (i, j, k, m = 1, ..., 3)
    [~, dG, d2G] = dChristMat(Cij, p, [0, 1, 1]);

    % Derivatives of the Christoffel equation F = det(G)
    for k = 1:3
        % f(k) = [d F]/[d p(k)]
        f(k) = dDetdX(G, dG(:,:,k));             
        for m = k:3
            % ff(k,m) = [d^2 F]/[d p(k) d p(m)]
            ff(k,m) = d2DetdXdY(G, dG(:,:,k), dG(:,:,m), d2G(:,:,k,m));  
        end;
    end;
    ff = symMat(ff);  % symmetrize matrix 'ff' 

    %% Construct the matrix Q = -f(3)^3*d^2 p(3)/[d p(k), d p(m)],  (k, m = 1, 2)
    for k = 1:2   
        for m = k:2
            Q(k,m) = ff(k,m)*f(3)^2 - ff(k,3)*f(m)*f(3) ...
                                    - ff(m,3)*f(k)*f(3) + ff(3,3)*f(k)*f(m);
        end;   
    end;
    Q = symMat(Q);  % (*) Matrix Q is singular, det(Q) = 0, at the shear-wave singularities

    %% Compute the NMO ellipse
    W = f(3)^2*dot(p, f)*inv(Q);  

    %% Expand the NMO ellipce to a cylinder
    for i = 1:2;
        q(i) = -f(i)/f(3);   
    end;
    Ucyl = [W(1,1),  W(1,2),  q(1)*W(1,1) + q(2)*W(1,2); ...
                 0,  W(2,2),  q(1)*W(1,2) + q(2)*W(2,2); ...
                 0,       0,  q(1)^2*W(1,1) + 2*q(1)*q(2)*W(1,2) + q(2)^2*W(2,2)];
    Ucyl = symMat(Ucyl);   
    
end;

end    % of the function
