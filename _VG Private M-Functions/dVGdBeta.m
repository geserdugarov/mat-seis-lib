%% *dVGdBeta*
% Derivatives of the phase velocity and group-velocity vector with respect to the phase angles

%%
% *Input:*

% Cij      - [6, 6] stiffness matrix
% n        - [3, 1] wavefront normal vector defined in terms of its spherical angles as
%            n = [sin(beta(1))*cos(beta(2)), sin(beta(1))*sin(beta(2)), cos(beta(1))]'  
% V        - [3, 1] phase velocities of three isonormal waves 
%            (typically sorted in function 'velPhaseU') 
% U        - [3, 3] polarizations (column vectors) of three isonormal waves 
%            (typically sorted in function 'velPhaseU') 
% waveType - [scalar] equal to 1 (for the P-wave), 2 (for the S1-wave or SV-wave), 
%            or 3 (for the S2-wave or SH-wave) that specifies the wave mode derivatives of 
%            whose velocities are to be computed  

%%
% *Output:*

% dVdb     - [1, 2] derivatives d V/d beta, where V is the phase velocity 
% dGdb     - [3, 2] derivatives d G/d beta, where G is the group-velocity vector

%%
% *Author:* Vladimir Grechka 2012 - 2014

%%
function [dVdb, dGdb] = dVGdBeta(Cij, n, V, U, waveType)
%% Settings  
[~, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

%% First and second-order derivatives of the elements of the Christoffel matrix with respect to beta
j1 = [1, 1, 2];   j2 = [1, 2, 2];   % indexes for the second-order derivatives of
                                    % the Christoffel matrix  Gamma = christoffel(Cij, n, n)    

% The first- and second-order derivatives of the wavefront normal n with respect to beta:
% dndb(:,j) = d n/d beta(j) and d2ndb2(:, j1+j2-1) = d^2 n/(d beta(j1) d beta(j2))
[dndb, spn, d2ndb2] = dVecdAngle(n);

dGamma = NaN(3, 3, 2);
for j = 1:2
    tmp = christoffel(Cij, n, dndb(:,j));           % d Gamma(:,:)/d beta(j) 
    dGamma(:,:,j) = tmp + tmp';
end;

d2Gamma = NaN(3, 3, 3);
for j = 1:3
    tmp1 = christoffel(Cij, n, d2ndb2(:,j));        % d^2 Gamma(:,:)/(d beta(j1) d beta(j2))
    tmp2 = christoffel(Cij, dndb(:,j1(j)), dndb(:,j2(j)));
    d2Gamma(:,:,j) = tmp1 + tmp1' + tmp2 + tmp2';
end;

%% First-order derivatives of the polarization vector with respect to beta

[~, dUdN] = dUdCijN(Cij, n, V, U, waveType, [0, 1, 0]); 
dU = dUdN*dndb;

%% First-order derivatives of the phase velocity
dVdb = NaN(2, 1);
for j = 1:2
    dV2     = U(:,waveType)'*dGamma(:,:,j)*U(:,waveType);                   % d V^2/d beta(j)
    dVdb(j) = dV2/(2*V(waveType));                                          % d V/d beta(j)
end;

%% Compute the derivatives d G/d beta if their output is required
d2Vdb2 = NaN(1, 3);
if nargout > 1    
    % Second-order derivatives of the phase velocity
    for j = 1:3
        d2V2 = 2*U(:,waveType)'*dGamma(:,:,j1(j))*dU(:,j2(j)) + ...         % d^2 V^2/(d beta(j1)*
                 U(:,waveType)'*d2Gamma(:,:,j)*U(:,waveType);               %          d beta(j2))
        % It is easy to verify that d2V2 is symmetric with respect to interchanging j1 and j2
        d2Vdb2(j) = (d2V2 - 2*dVdb(j1(j))*dVdb(j2(j)))/(2*V(waveType));     % d^2 V/(d beta(j1)*
    end;                                                                    %        d beta(j2))
        
    % Derivatives of the group velocity d G(:)/d beta(j) derived in Grechka and Obolentseva 
    % (1989, p 19, eq 33) 
    dGdb(:,1) = dndb(:,1)*(V(waveType) + d2Vdb2(1)) + ...
                dndb(:,2)*(d2Vdb2(2) - n(3)*dVdb(2)/spn)/spn^2;
    dGdb(:,2) = n*dVdb(2) + dndb(:,1)*d2Vdb2(2) + dndb(:,2)*(V(waveType) + d2Vdb2(3)/spn^2) + ...
                d2ndb2(:,2)*dVdb(1) + d2ndb2(:,3)*dVdb(2)/spn^2;
end;
    
end    % of the function