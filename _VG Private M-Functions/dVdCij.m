%% *dVdCij*
% Derivatives of the phase velocity with respect to the stiffness coefficients derived in
% Grechka and Duchkov (2011)

%%
% *Input:*

% n     - [3, 1] wavefront normal vector 
% V     - [scalar] phase velocity
% U     - [3, 1] polarization vector

%%
% *Output:*

% dVmat - [6, 6]  matrix of the Frechet derivatives dV/dC(i,j) 
% dVvec - [1, 21] vector of the Frechet derivatives dV/dC(i,j) 

%%
% *Author:* Vladimir Grechka 2012 

%%
function [dVmat, dVvec] = dVdCij(n, V, U) 
%% Settings  
[~, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

[~, ~, ~, t2m] = indexesVMT(4);      % get the indexes

%% Implementation of Appendix A in Grechka and Duchkov (2011)
dV2 = zeros(6, 6);
for i = 1:3   
    for j = 1:3   
        for k = 1:3   
            for l = 1:3
                im = t2m(i,j);   jm = t2m(k,l); 
                f = 1;
                if im ~= jm
                    f = 2;           % account for the symmetry d/dCij(im,jm) = d/dCij(jm,im)
                end; 
                dV2(im,jm) = dV2(im,jm) + f*U(i)*n(j)*n(k)*U(l);   % eq A-5 in Grechka and Duchkov (2011)
            end;   
        end;   
    end;   
end;

dVmat = dV2/(2*V);
dVvec = mat2vec(dVmat);

end    % of the function