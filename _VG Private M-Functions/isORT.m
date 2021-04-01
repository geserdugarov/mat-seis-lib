%% *isORT*
% Check whether input stiffness matrix Cij is orthorhombic (ORT) within given tolerance

%%
% *Input:*

% Cij  - [6, 6] stiffness matrix
% tol  - [scalar] tolerance of proximity to orthotropy 

%%
% *Output:*

% flag - [scalar] equal to 1 if Cij is ORT within given tolerance and to 0 otherwise
% R    - [3, 3] rotation matrix that transforms Cij into its ORT form

%%
% *Author:* Vladimir Grechka 2014

%%
function [flag, R] = isORT(Cij, tol)
%% Settings and defaults
% [thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
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
[Dvec, Dval] = eig(symMat(Dij));                    % eigenvectors and eigenvalues of 'Dij'

% Check whether 'Dij' has double eigenvalues
dDval = [abs(Dval(2,2) - Dval(3,3)), abs(Dval(1,1) - Dval(3,3)), abs(Dval(1,1) - Dval(2,2))]; 
[Dvalsort, isort] = sort(dDval, 'descend'); 
Dvec = Dvec(:, isort);   % place possibly double eigenvectors in the first two columns of 'Dvec'
if Dvalsort(3) < tol
    doubleIndexDij = 1;                             % 'Dij' has double eigenvalues
else
    doubleIndexDij = 0;                             % 'Dij' has distinct eigenvalues
end;

%% Voigt stiffness tensor (Helbig, 1994, p 111; Cowin, 1989)
Vij = zeros(3, 3);
for i = 1:3   
   for j = i:3
       Vij(i,j) = Cij(t2m(i,1), t2m(j,1)) + Cij(t2m(i,2), t2m(j,2)) + Cij(t2m(i,3), t2m(j,3));
   end;   
end;
[Vvec, Vval] = eig(symMat(Vij));                    % eigenvectors and eigenvalues of 'Vij'

% Check whether 'Vij' has double eigenvalues
dVval = [abs(Vval(2,2) - Vval(3,3)), abs(Vval(1,1) - Vval(3,3)), abs(Vval(1,1) - Vval(2,2))];
[Vvalsort, jsort] = sort(dVval, 'descend');
Vvec = Vvec(:, jsort);   % place possibly double eigenvectors in the first two columns of 'Vvec'
if Vvalsort(3) < tol
    doubleIndexVij = 1;                             % 'Vij' has double eigenvalues
else
    doubleIndexVij = 0;                             % 'Vij' has distinct eigenvalues
end;

% %% Hydrostatic pressure compliance tensor (Helbig, 1994, p 111; Cowin, 1989)
% Sij = inv(Cij);
% for i = 1:3   
%     for j = i:3
%         Vij(i,j) = Sij(t2m(i,j),1) + Sij(t2m(i,j),2) + Sij(t2m(i,j),3);
%     end;   
% end;
% Vij = symMat(Vij).*[1, 0.5, 0.5; 0.5, 1, 0.5; 0.5, 0.5, 1];
% [Vvec, ~] = eig(Vij);                               % eigenvectors and eigenvalues of Vij

%% Change the sign of 'Dvec' if necessary
if dot(Dvec(:,3), Vvec(:,3)) < 0
    Dvec = -Dvec;
end;

%% Rotate the eigenvectors if necessary 
if doubleIndexDij == 1  ||  doubleIndexVij == 1
    % The case of double eigenvalues of either 'Dij' or 'Vij':
    % Attempt to rotate eigenvectors in the first two columns of 'Dvec' to
    % align them with the eigenvectors in the first two columns of 'Vvec'
    rotMat1  = Dvec'*Vvec;
    rotAngle = atan2(rotMat1(2,1), rotMat1(1,1));       % rotation angle around 'Dvec(:,3)'
    rotMat2  = rotateAroundAxis(Dvec(:,3), rotAngle, 'rad');
    Dvec     = rotMat2*Dvec;                            % the rotated 'Dvec'
else
    % The case of distinct eigenvalues of both 'Dij' and 'Vij':
    % If 'Dij' and 'Vij' have all distinct eigenvalues, the equality of eigenvectors 
    % 'Dvec' and 'Vvec' is a necessary criterion for orthotropy -- no action is needed  
end;    
    
%% Compare the eigenvectors of 'Dij' and 'Vij'
count = 0;    
for i = 1:3
    for j = 1:3
%        norm(cross(Dvec(:,i), Vvec(:,j)))
        if norm(cross(Dvec(:,i), Vvec(:,j))) < tol
            count = count + 1;
        end;
    end;
end;

if count == 3
    R = Dvec;                                        
    CijRot = bond(Cij, R);                           % rotation to ORT if 'Cij' is a tilted ORT
    
    % Check whether 'CijRot' is, indeed, an ORT matrix
    check = sum(abs(CijRot(1:3, 4))) + ...          % off-diagonal zeros
            sum(abs(CijRot(1:4, 5))) + ...
            sum(abs(CijRot(1:5, 6)));   
    if check < tol                                     
        flag = 1;                                   % 'Cij' has ORT symmetry
    else
        flag = 0;                                   % The symmetry of 'Cij' is lower than ORT
        R = eye(3, 3);
    end;

else                                                 
    flag = 0;                                       % The symmetry of 'Cij' is lower than ORT
    R = NaN(3, 3);
end;

end    % of the function
