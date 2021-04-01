%% *isMNC*
% Check whether input stiffness matrix Cij is monoclinic (MNC) within given tolerance

%%
% *Input:*

% Cij  - [6, 6] stiffness matrix
% tol  - [scalar] tolerance of proximity to monoclinic symmetry 

%%
% *Output:*

% flag - [scalar] equal to 1 if Cij is MNC within given tolerance and to 0 otherwise
% R    - [3, 3] rotation matrix that transforms Cij into its MNC form

%%
% *Author:* Vladimir Grechka 2014

%%
function [flag, R] = isMNC(Cij, tol)
%% Settings and defaults
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
% narginTrue = nargin(thisFileName);   
% narginchk(narginTrue, narginTrue);

if nargin < 2;   tol = 1.e-6;   end;                % default for tol

%% Dilatational stiffness tensor (Helbig, 1994, p 111; Cowin, 1989)
Dij = zeros(3, 3);
[~, ~, ~, t2m] = indexesVMT(4);                     % get the indexes
for i = 1:3   
    for j = i:3
        Dij(i,j) = Cij(t2m(i,j),1) + Cij(t2m(i,j),2) + Cij(t2m(i,j),3);
    end;   
end;
[Dvec, ~] = eig(symMat(Dij));                       % eigenvectors and eigenvalues of 'Dij'

% %% Voigt stiffness tensor (Helbig, 1994, p 111; Cowin, 1989)
% Vij = zeros(3, 3);
% for i = 1:3   
%    for j = i:3
%        Vij(i,j) = Cij(t2m(i,1), t2m(j,1)) + Cij(t2m(i,2), t2m(j,2)) + Cij(t2m(i,3), t2m(j,3));
%    end;   
% end;
% [Vvec, ~] = eig(symMat(Vij))                       % eigenvectors and eigenvalues of 'Vij'

%% Hydrostatic pressure compliance tensor (Helbig, 1994, p 111; Cowin, 1989)
Vij = zeros(3, 3);
Sij = inv(Cij);
for i = 1:3   
    for j = i:3
        Vij(i,j) = Sij(t2m(i,j),1) + Sij(t2m(i,j),2) + Sij(t2m(i,j),3);
    end;   
end;
Vij = symMat(Vij).*[1, 0.5, 0.5; 0.5, 1, 0.5; 0.5, 0.5, 1];
[Vvec, ~] = eig(Vij);                               % eigenvectors and eigenvalues of 'Vij'

%% Compare the eigenvectors of Dij and Vij
% Identity of at least one eigenvector of Dij and Vij is the criterion for monoclinic symmetry
count = 0;  ii = [];
for i = 1:3
    for j = 1:3
        if norm(cross(Dvec(:,i), Vvec(:,j))) < tol
            count = count + 1;
            ii = i; % 'ii' is the axis normal to the symmetry plane that will become x3 if 
                    % 'Dij' and 'Vij' have a single common egenvector, that is, final 'count' = 1 
        end;           
    end;
end;

if count >= 1
    flag = 1;                                       % Cij has MNC symmetry or higher
    R = Dvec; 
    if count == 1
        if ii == 1                                  % If ii ~= 3, permute the axes to make
            R = [Dvec(:,2), Dvec(:, 3), Dvec(:,1)]; % the symmetry plane horizontal     
        elseif ii == 2                                  
            R = [Dvec(:,3), Dvec(:, 1), Dvec(:,2)];
        end;
        % Add rotation in the [x1, x2]-plane  
        CijR = bond(Cij, R);
        theta = (1/2)*atan2( 2*CijR(4,5), CijR(5,5) - CijR(4,4));  % Eq A3 in Grechka et al. (2000)
        Rx1x2 = [cos(theta), sin(theta), 0;  -sin(theta), cos(theta), 0;  0, 0, 1];
        R = R*Rx1x2';                               % This rotation can be verified 
    else                                            %   with 'CijRot = bond(Cij, R)'
        fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
        display('>>> WARNING: Two chosen rank two contractions of the stiffness tensor');
        display(Cij);
        display('    have more than one common eigenvector ');
        fprintf('>>> Rotation to a horizontal symmetry plane is expected to be ambigous \n');
    end;
else
    flag = 0;                                       % The symmetry of Cij is lower than MNC
    R = eye(3, 3);
end;

end    % of the function
