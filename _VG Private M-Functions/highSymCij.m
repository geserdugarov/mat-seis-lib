%% *highSymCij*
% Find high-symmetry approximations (isotropic, VTI, and orthorhombic) of a given stiffness matrix

%%
% *Input:*

% Cij  - [6, 6] stiffness matrix  

%%
% *Output:*

% Ciso - [6, 6] isotropic approximation of Cij
% Cvti - [6, 6] VTI approximation of Cij
% Cort - [6, 6] orthorhombic approximation of Cij

%%
% *Author:* Vladimir Grechka 2012-2014

%% 
% *Known issues:* 
%
% * Matrix Cij should be brought as close as possible to a crystallographic coordinate frame
%   before applying 'highSymCij'

%%
function [Ciso, Cvti, Cort] = highSymCij(Cij)
%% Settings 
[~, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

%% Isotropic approximation (Fedorov, 1968) 
sum0 = Cij(1,1) + Cij(2,2) + Cij(3,3);
sum1 = sum0 + 2*(Cij(1,2) + Cij(1,3) + Cij(2,3));
sum2 = sum0 + 2*(Cij(4,4) + Cij(5,5) + Cij(6,6));
lambda = (2*sum1 - sum2)/15;   mu = (3*sum2 - sum1)/30;  
Ciso = thomsen2cij([sqrt(lambda + 2*mu), sqrt(mu), 0, 0, 0]); 

%% VTI approximation (Dellinger, 2005, equation 3)
if nargout > 1
    Cvti = zeros(6,6);
    Cvti(1,1) = (3*Cij(1,1) + 2*Cij(1,2) + 3*Cij(2,2) + 4*Cij(6,6))/8;
    Cvti(2,2) = Cvti(1,1);
    Cvti(1,2) = (Cij(1,1) + 6*Cij(1,2) + Cij(2,2) - 4*Cij(6,6))/8;
    Cvti(1,3) = (Cij(1,3) + Cij(2,3))/2;
    Cvti(2,3) = Cvti(1,3);
    Cvti(3,3) = Cij(3,3);
    Cvti(4,4) = (Cij(4,4) + Cij(5,5))/2;
    Cvti(5,5) = Cvti(4,4);
    Cvti(6,6) = (Cij(1,1) - 2*Cij(1,2) + Cij(2,2) + 4*Cij(6,6))/8;
    Cvti = symMat(Cvti);
end;

%% Orthorhombic approximation obtained from the requirement norm(Cij - Cort) = min
if nargout > 2
    Cort = zeros(6,6);
    for i=1:3
        Cort(i+3,i+3) = Cij(i+3,i+3);
        for j=i:3
            Cort(i,j) = Cij(i,j);
        end;
    end;
    Cort = symMat(Cort);
end;

end    % of the function