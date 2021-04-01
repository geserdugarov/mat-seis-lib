%% *dGdCij*
% Derivatives of the group velocity with respect to the stiffness coefficients derived in
% Grechka and Duchkov (2011)

%%
% *Input:*

% Cij   - [6, 6] stiffness matrix
% p     - [3, 1] slowness vector 
% U     - [3, 1] polarization vector
% g     - [3, 1] group-velocity vector

%%
% *Output:*

% dGmat - [6, 6]  matrix of the Frechet derivatives dV/dC(i,j) 
% dGvec - [1, 21] vector of the Frechet derivatives dV/dC(i,j) 

%%
% *Author:* Vladimir Grechka 2012 - 2014

%%
function [dGmat, dGvec] = dGdCij(Cij, p, U, g)
%% Settings  
[~, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

tol = 1.e-12;

%% Derivatives of the phase velocity
V = 1/norm(p);
n = V*p;
[dVmat, dVvec] = dVdCij(n, V, U);

if isISO(Cij, tol) == 1
    %% Isotropy
    dGmat = dVmat;                         
    dGvec = dVvec;
    
else
    %% Anisotropy
    G = christoffel(Cij, p, p) - eye(3);        % construct the Christoffel matrix G(i,j)
    [~, dG, ~] = dChristMat(Cij, p, [0, 1, 0]); % differentiate its elements d G(i,j)/d p(m)
            
    Fp = zeros(3, 1);
    for m = 1:3
        Fp(m) = dDetdX(G, dG(:,:,m));      % differentiate the determinant dF = d det(G)/dp(m)
    end;
        
    if norm(Fp) > tol
        % Calculate dg/dCij (Grechka and Duchkov, 2011, eq A-14)
        dFdC = dChristdCij(Cij, p, V);      % differentiate det of Christoffel matrix dF/dC(i,j)
        dGmat = dot(g, g)*dFdC/sqrt(dot(Fp, Fp));
        dGvec = mat2vec(dGmat);
        
        % Make sure that the signs of dG and dV coincide
        iv = find(abs(dVvec) == max(abs(dVvec)));
        sgn = sign(dVvec(iv(1))*dGvec(iv(1)));
        dGmat = sgn*dGmat;
        dGvec = sgn*dGvec;    
    else
        % Equality |Fp| = 0 indicates the shear-wave singularity at which dG/dCij blows up 
        dGmat = dVmat;     % --> approximate dG/dCij with dV/dCij
        dGvec = dVvec;
    end;
      
end;

end    % of the function

%%