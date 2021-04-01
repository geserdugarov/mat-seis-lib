%% *dUdCijN*
% Derivatives of components of the polarization vector with respect to the stiffness coefficients 
% and components of the wavefront normal vector

%%
% *Input:*

% Cij      - [6, 6] stiffness matrix in Voigt notation
% n        - [3, 1] wavefront normal vector
% V        - [3, 1] phase velocities of three isonormal waves 
%            (typically sorted in function 'velPhaseU') 
% U        - [3, 3] polarizations (column vectors) of three isonormal waves 
%            (typically sorted in function 'velPhaseU') 
% waveType - [scalar] equal to 1 (for the P-wave), 2 (for the S1-wave or SV-wave), 
%            or 3 (for the S2-wave or SH-wave) that specifies the wave mode whose polarization
%            derivatives are to be computed  
% flag     - [1, 2] array whose values [1, 0], [0, 1], or [1, 1] determine whether either
%            output array (dGdC or dGdN) or both array (dGdC and dGdN) are to be computed

%%
% *Output:*

% dUdC     - [3, 21] array of the derivatives dU(i,waveType)/dCij
% dUdN     - [3, 3]  array of the derivatives dU(i,waveType)/dn(k)

%%
% *Author:* Vladimir Grechka 2012 2013

%%
function [dUdC, dUdN] = dUdCijN(Cij, n, V, U, waveType, flag) 
%% Settings and checks 
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

% Check whether the value of variable flag is legitimate  
if isempty(flag) == 1  ||  size(flag, 1) ~= 1  ||  size(flag, 2) ~= 3
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Incorrectly specified variable ''flag'' \n');
    display(flag)
      error('>>> STOP');   
end;

% Check whether variable waveType contains legitimate values 1, 2, or 3
if isempty(find((waveType ~= 1  &  waveType ~= 2  &  waveType ~= 3), 1)) == 0
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Array ''waveType'' contains elelemt(s) different from 1, 2, or 3 \n');
    display(waveType)    
    fprintf('>>> PAUSE \n');    pause;
end;   

% Check whether the wavefront normal is a unit vector
tol = 1.e-6; %tol = 1.e-14;
if abs(norm(n) - 1) > tol
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Wavefront normal ''n'' is not a unit vector within specified tolerance \n');
    fprintf('>>> abs(norm(n) - 1) = %g \n \n', abs(norm(n) - 1));
      error('>>> STOP');   
end;

% Check whether the polarization matrix is unitary
if isUnitary(U, tol) == 0
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> The polarization matrix ''U'' is not unitary within specified tolerance \n');
    fprintf('>>> U*U'' - I = \n');
    display(U*U' - eye(3));
      error('>>> STOP');   
end;

%% Detect the shear-wave singularity
dV2 = V(waveType)^2 - V.^2;
dV2sorted = sort(abs(dV2), 'ascend');

if dV2sorted(2)/V(waveType)^2 < tol  ||  dV2sorted(3)/V(waveType)^2 < tol
    dV2(2) =  sqrt(tol);
    dV2(3) = -sqrt(tol);
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Shear-wave singularity: V = [%g, %g, %g] at n = [%g, %g, %g] \n', [V, n']);
    fprintf('>>> Division by zero in computing the polarization derivatives is avoided \n');
    fprintf('    by taking the default values dV2(2) = %g and dV2(3) = %g \n', dV2(2), dV2(3));
end;   

%% Derivatives of the Christoffel matrix
[dGdC, dGdN, ~] = dChristMat(Cij, n, flag);
[v2m, ~, ~, ~] = indexesVMT(1);

%% dU/dC
dUdC = zeros(3,21);
if flag(1,1) == 1
    for i = 1:21
        dG = dGdC(:,:,v2m(1,i),v2m(2,i));
        for iwave = 1:3
            if iwave ~= waveType
                dUdC(:,i) = dUdC(:,i) + (U(:,iwave)'*dG*U(:,waveType)/dV2(iwave))*U(:,iwave);
            end;
        end;
    end;
end;

%% dU/dn
dUdN = zeros(3,3);
if flag(1,2) == 1
    for i = 1:3
        dG = dGdN(:,:,i);
        for iwave = 1:3
            if iwave ~= waveType
                dUdN(:,i) = dUdN(:,i) + (U(:,iwave)'*dG*U(:,waveType)/dV2(iwave))*U(:,iwave);
            end;
        end;
    end;
end;

end    % of the function