%% *nmoEll2Cyl*
% Reconstruct the NMO cylinder from its cross-section by a plane (Grechka and Tsvankin, 2002;
% Appendixes B, D)

%%
% *Input:*
% W0   - [2, 2] matrix expressing the NMO ellipse, which is a cross-section of the NMO cylinder 
%        by a plane
% nrm  - [3, 1] vector normal to the cross-section plane
% q    - [2, 1] vector of the slowness derivatives dp(3)/dp(i), (i = 1, 2) computed by 
%        function 'nmoCylInt'

%%
% *Output:*
% Ucyl - [3, 3] symmetric matrix representing the NMO cylinder

%%
% *Author:* Vladimir Grechka 1998, 2014

%%
% *Comments:*
%
% * Original Matlab version, called 'ell2cyl', is a part of the 'ART' package freely  
%   distributed by the CWP, CSM

%%
function [Ucyl] = nmoEll2Cyl(W0, nrm, q)
%% Settings  
[~, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

Dmat = zeros(2,2,2,2);

%% Construct vectors normal to 'nrm' that lie in the plane of the NMO ellipse
[azm, pol90, ~] = cart2sph(nrm(1), nrm(2), nrm(3));
pol = pi/2 - pol90;
b(:,1) = [cos(pol)*cos(azm), cos(pol)*sin(azm), -sin(pol)]';
b(:,2) = [-sin(azm), cos(azm), 0]';

%%
% Build the 3*3 matrix to be inverted to find the NMO ellipse 'W', which is the cross-section
% of the sought NMO cylinder 'Ucyl' by the horizontal plane (equations D7 and D9 in 
% Grechka and Tsvankin, 2002)
for i = 1:2
    for j = 1:2
        for k = 1:3
            for m = 1:3
                Btmp = (b(k,i)*b(m,j) + b(k,j)*b(m,i))/2;
                if k == 1  &&  m == 1
                    Dmat(i,j,1,1) = Dmat(i,j,1,1) + Btmp; 
                elseif k == 1  &&  m == 2
                    Dmat(i,j,1,2) = Dmat(i,j,1,2) + Btmp;
                elseif k == 1  &&  m == 3
                    Dmat(i,j,1,1) = Dmat(i,j,1,1) + q(1)*Btmp;
                    Dmat(i,j,1,2) = Dmat(i,j,1,2) + q(2)*Btmp;
                elseif k == 2  &&  m == 1
                    Dmat(i,j,2,1) = Dmat(i,j,2,1) + Btmp;
                elseif k == 2  &&  m == 2
                    Dmat(i,j,2,2) = Dmat(i,j,2,2) + Btmp;
                elseif k == 2  &&  m == 3
                    Dmat(i,j,1,2) = Dmat(i,j,1,2) + q(1)*Btmp;
                    Dmat(i,j,2,2) = Dmat(i,j,2,2) + q(2)*Btmp;
                elseif k == 3  &&  m == 1
                    Dmat(i,j,1,1) = Dmat(i,j,1,1) + q(1)*Btmp;
                    Dmat(i,j,2,1) = Dmat(i,j,2,1) + q(2)*Btmp;
                elseif k == 3  &&  m == 2
                    Dmat(i,j,2,1) = Dmat(i,j,2,1) + q(1)*Btmp;
                    Dmat(i,j,2,2) = Dmat(i,j,2,2) + q(2)*Btmp;
                else    %  k == 3  &&  m == 3
                    Dmat(i,j,1,1) = Dmat(i,j,1,1) + q(1)^2*Btmp;
                    Dmat(i,j,1,2) = Dmat(i,j,1,2) + q(1)*q(2)*Btmp;
                    Dmat(i,j,2,1) = Dmat(i,j,2,1) + q(1)*q(2)*Btmp;
                    Dmat(i,j,2,2) = Dmat(i,j,2,2) + q(2)^2*Btmp;
                end; 
            end;
        end;
    end;
end;

%% Equation D11 in Grechka and Tsvankin (2002)
Emat = [Dmat(1,1,1,1), 2*Dmat(1,1,1,2), Dmat(1,1,2,2); ...
     Dmat(1,2,1,1), 2*Dmat(1,2,1,2), Dmat(1,2,2,2); ...
     Dmat(2,2,1,1), 2*Dmat(2,2,1,2), Dmat(2,2,2,2)];   
Wnrm = [W0(1,1), W0(1,2), W0(2,2)]';
W = Emat\Wnrm;

%% Reconstruct the NMO cylinder
Ucyl = [W(1),  W(2),  q(1)*W(1) + q(2)*W(2); ...
           0,  W(3),  q(1)*W(2) + q(2)*W(3); ...
           0,     0,  q(1)^2*W(1) + 2*q(1)*q(2)*W(2) + q(2)^2*W(3)];
Ucyl = symMat(Ucyl);   

end    % of the function
