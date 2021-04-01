%% *tsvankin2cij*
% Calculate stiffness matrix in an orthorhombic medium from its Tsvankin's parameters 

%%
% *Input:*

% Ani - [1, 9] vector of Tsvankin's anisotropy parameters arranged as
%       [Vp0, Vs0, epsilon1, epsilon2, delta1, delta2, delta3, gamma1, gamma2]

%%
% *Output:*

% Cij - [6, 6] ORT stiffness matrix

%%
% *Author:* Vladimir Grechka 1996 2012 - 2014

%%
function [Cij] = tsvankin2cij(Ani)
%% Settings and defaults
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

%% Calculate the stiffnesses
c33 = Ani(1)^2;   c22 = c33*(2*Ani(3) + 1);   c11 = c33*(2*Ani(4) + 1);      
c55 = Ani(2)^2;   c66 = c55*(2*Ani(8) + 1);   c44 = c66/(2*Ani(9) + 1);

if c33 < c55   
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Untypical parameters encountered: c33 (= %g) < c55 (= %g) \n', c33, c55);
end;

cc = (c33 - c44)*(2*Ani(5)*c33 + c33 - c44);
if cc > 0   
    if c33 > c44
        c23 = sqrt(cc) - c44;
    else
        c23 = -sqrt(cc) - c44;
        fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
        fprintf('>>> Untypical parameters encountered: c33 (= %g) < c44 (= %g) \n', c33, c44);
    end;
else
    c23 = 0;
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Coefficient c23 is set to zero to avoid a complex-valued stiffness matrix \n');
end;

cc = (c33 - c55)*(2*Ani(6)*c33 + c33 - c55);
if cc > 0   
    if c33 > c55
        c13 = sqrt(cc) - c55;
    else
        c13 = -sqrt(cc) - c55;
        fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
        fprintf('>>> Untypical parameters encountered: c33 (= %g) < c55 (= %g) \n', c33, c55);
    end;
else 
    c13 = 0;
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Coefficient c13 is set to zero to avoid a complex-valued stiffness matrix \n');
end;

cc = (c11 - c66)*(2*Ani(7)*c11 + c11 - c66);
if cc > 0   
    if c11 > c66
        c12 = sqrt(cc) - c66;
    else
        c12 = -sqrt(cc) - c66;
        fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
        fprintf('>>> Untypical parameters encountered: c11 (= %g) < c66 (= %g) \n', c11, c66);  
    end;
else
    c12 = 0;
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Coefficient c12 is set to zero to avoid a complex-valued stiffness matrix \n');
end;

%% Construct the stiffness matrix
Cij = zeros(6, 6);
Cij(1,1) = c11;   Cij(1,2) = c12;   Cij(1,3) = c13;
                  Cij(2,2) = c22;   Cij(2,3) = c23;   Cij(3,3) = c33;
Cij(4,4) = c44;   Cij(5,5) = c55;   Cij(6,6) = c66;
Cij = symMat(Cij);   % fill the symmetric part

end   % of the function