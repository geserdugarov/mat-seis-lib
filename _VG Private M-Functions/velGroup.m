%% *velGroup*
% Compute the group velocity using the algorithm of Grechka and Duchkov (2011) 

%%
% *Input:*

% Cij      - [6, 6] stiffness matrix in Voigt notation
% n        - [3, 1] wavefront normal 
% V        - [3, 1] phase velocities of three isonormal waves (typically sorted in 
%            function 'velPhaseU') 
% U        - [3, 3] polarizations (column vectors) of three isonormal waves (also sorted in 
%            'velPhaseU')  
% waveType - vector containing a sequence of numbers 1 (for the P-wave), 2 (for the S1-wave or
%            SV-wave), or 3 (for the S2-wave or SH-wave) that specify the wave modes for which
%            the group velocities are to be computed  

%%
% *Output:*

% g        - [3, 3] matrix of the group velocities of isonormal waves defined by 'waveType'

%%
% *Author:* Vladimir Grechka 2012 - 2014

%%
function [g] = velGroup(Cij, n, V, U, waveType)
%% Settings 
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

g = NaN(3,3);
[~, ~, ~, t2m] = indexesVMT(4);

%% Various checks
% Check whether 'waveType' contains ligitimate values 1, 2, or 3
if isempty(find((waveType ~= 1  &  waveType ~= 2  &  waveType ~= 3), 1)) == 0
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Array ''waveType'' contains elelemt(s) different from 1, 2, or 3 \n');
    display(waveType)    
    fprintf('>>> PAUSE \n');    pause;
end;   

% Check whether the wavefront normal is a unit vector
tol = 1.e-8;
if abs(norm(n) - 1) > tol
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Wavefront normal ''n'' is not a unit vector within specified tolerance \n');
    fprintf('>>> abs(norm(n) - 1) = %g \n', abs(norm(n) - 1));
      error('>>> STOP');   
end;

%% Compute the group velocity
if isISO(Cij) == 1
    for iw = 1:3
        g(:,iw) = V(iw)*n;      % isotropy
    end;
else
    % Check whether the polarization matrix is unitary
    if isUnitary(U, tol) == 0
        fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
        fprintf('>>> The polarization matrix ''U'' is not unitary within specified tolerance \n');
        fprintf('>>> U*U'' - I = \n');
        display(U*U' - eye(3))
    end;
    
    % Rotate the wavefront normal vector and the stiffness matrix to the coordinate frame, 
    % in which U = eye(3) (Grechka and Duchkov, 2011)
    nRot = U'*n;
    CijRot = bond(Cij, U);

    for i = 1:length(waveType)
        iwave = waveType(i);
        pRot  = nRot/V(iwave);                       % the slowness vector in the rotated frame
        CijQ  = CijRot(t2m(iwave,:), t2m(iwave,:));  % eq 19 in Grechka and Duchkov (2011)
        gRot  = CijQ*pRot;                           % the group velocity in the rotated frame        
        g(:,iwave) = U*gRot;                         % the group velocity in the original frame
    end;
end;

end    % of the function