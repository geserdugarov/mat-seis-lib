%% *longitudinalDir*
% Compute longitudinal directions, defined as Up = n, for a given stiffness matrix

%%
% *Input:*

% Cij - [6, 6] stiffness matrix

%%
% *Output:*

% nL - array of the unit longitudinal normals

%%
% *Author:* Vladimir Grechka 20134

%%
function [nL] = longitudinalDir(Cij)
%% Settings  
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

syms n1 n2 n3;
tol1 = 1.e-12;  tol2 = 1.e-6;

%% Check whether the model is isotropic
if isISO(Cij) == 1
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    display('>>> Isotropic stiffness tensor Cij = ');
       disp(Cij);
    fprintf('>>> Any wavefront normal is a singularity in isotropic media -- PAUSE \n');  pause;
    nL = NaN(1, 3);  
    return;
end;

%% TI and low-symmetry media
[flag, R] = isTI(Cij, tol1);
if flag == 1
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    display('>>> The stiffness tensor');
       disp(Cij);
    fprintf('>>> is TI within the tolerance %g -- Look for solution in the [n1, n3]-plane  \n', tol1); 
    C1 = bond(Cij, R);
    
    % The Christoffel matrix G = n.Cij.n in the [n1, n3]-plane for the P- and SV-waves
    G11 = C1(1,1)*n1^2 + C1(5,5)*n3^2;
    G33 = C1(5,5)*n1^2 + C1(3,3)*n3^2;
    G13 = (C1(1,3) + C1(5,5))*n1*n3;
    
    % Components of the dot product G.n 
    d1 = G11*n1 + G13*n3; 
    d3 = G13*n1 + G33*n3; 

    %% The cross product [n, G.n] = 0 (Fedorov, eq 18.2)
    digits(30);
    eq = n1*d3 - n3*d1;
    F1 = simplify(eq/(n1*n3));  % remove known roots [1, 0] and [0, 1] to speed up computations
    F0 = n1^2 + n3^2 - 1;

    %% Solve the system
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    display('>>> Computation of the longitudinal normals normally takes a few seconds...');
%    tic
    [s1, s3] = vpasolve([F0 == 0, F1 == 0], [n1, n3], [-1, 1; -1, 1]);
%    toc
    
    %% Add the known roots to the solution
    nTI(1,:) = [0, 0, 1];
    nTI(2,:) = [1, 0, 0];
    count = 2;    
    %% Select solutions that yield the equal n and U
    for isol = 1:length(s1)
        tmp = double([s1(isol), 0, s3(isol)]');  % convert symbolic variables to numeric form
        n = tmp/norm(tmp);
        [~, U] = velPhaseU(Cij, n, tol2);
        if abs( abs(dot(n, U(:,1))) - 1 ) < tol2  
            count = count + 1;
            nTI(count, :) = n';
        end;    
    end;

    for i = 1:size(nTI, 1)
        nL(i,:) = (R*nTI(i,:)')';
    end;
    
else
    %% Symmetries lower than TI
    % The Christoffel matrix G = n.Cij.n
    G11 = Cij(1,1)*n1^2 + Cij(6,6)*n2^2 + Cij(5,5)*n3^2 + 2*Cij(1,6)*n1*n2 + ...
                                                          2*Cij(1,5)*n1*n3 + 2*Cij(5,6)*n2*n3;
    G22 = Cij(6,6)*n1^2 + Cij(2,2)*n2^2 + Cij(4,4)*n3^2 + 2*Cij(2,6)*n1*n2 + ...
                                                          2*Cij(4,6)*n1*n3 + 2*Cij(2,4)*n2*n3;
    G33 = Cij(5,5)*n1^2 + Cij(4,4)*n2^2 + Cij(3,3)*n3^2 + 2*Cij(4,5)*n1*n2 + ...
                                                          2*Cij(3,5)*n1*n3 + 2*Cij(3,4)*n2*n3;
    G12 = Cij(1,6)*n1^2 + Cij(2,6)*n2^2 + Cij(4,5)*n3^2 + ...
         (Cij(1,2) + Cij(6,6))*n1*n2 + (Cij(1,4) + Cij(5,6))*n1*n3 + (Cij(4,6) + Cij(2,5))*n2*n3;
    G13 = Cij(1,5)*n1^2 + Cij(4,6)*n2^2 + Cij(3,5)*n3^2 + ...
         (Cij(1,4) + Cij(5,6))*n1*n2 + (Cij(1,3) + Cij(5,5))*n1*n3 + (Cij(4,5) + Cij(3,6))*n2*n3;
    G23 = Cij(5,6)*n1^2 + Cij(2,4)*n2^2 + Cij(3,4)*n3^2 + ...
         (Cij(4,6) + Cij(2,5))*n1*n2 + (Cij(4,5) + Cij(3,6))*n1*n3 + (Cij(2,3) + Cij(4,4))*n2*n3;
    
    % Components of the dot product G.n 
    d1 = G11*n1 + G12*n2 + G13*n3; 
    d2 = G12*n1 + G22*n2 + G23*n3; 
    d3 = G13*n1 + G23*n2 + G33*n3; 

    %% Components of the cross product [n, G.n] = 0 (Fedorov, eq 18.2)
    digits(30);
    F1 = n2*d3 - n3*d2;
    F2 = n1*d3 - n3*d1;
    F3 = n1*d2 - n2*d1;
    F0 = n1^2 + n2^2 + n3^2 - 1;

    %% Solve the system
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    display('>>> Computation of the longitudinal normals normally takes a few seconds...');
%    tic
    [s1, s2, s3] = vpasolve([F0 == 0, F1 == 0, F2 == 0, F3 == 0], [n1, n2, n3], ...
        [-1, 1; -1, 1; -1, 1]);
%    toc
        
    %% Select solutions that yield the equal n and U
    count = 0;
    for isol = 1:length(s1)
        tmp = double([s1(isol), s2(isol), s3(isol)]');  % convert symbolic variables to numeric form
        n = tmp/norm(tmp);
        [~, U] = velPhaseU(Cij, n, tol2);
        if abs( abs(dot(n, U(:,1))) - 1 ) < tol2  
            count = count + 1;
            nL(count, :) = n';
        end;    
    end;
    
end;

end  % of the function
