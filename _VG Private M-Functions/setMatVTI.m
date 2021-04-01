%% *setMatVTI*
% Compute matrices relating the slowness components and tangents of the group and phase angles
% of the P- and SV-waves propagating in VTI media  

%% 
% *Input:*  

% c11, c13, c33, c55 - stiffness coefficients governing P- and SV-wave propagation in VTI media

%% 
% *Output:*

% F - [3, 3] matrix representing the Christoffel equation in the form
%     [1, p1^2, p1^4] * F * [1, p3^2, p3^4]' = 0, 
%     where p = [p1, 0, p3] is the slowness vector
% M - [7, 3] matrix relating Psi = tan(group angle) to Theta = tan(phase angle) as
%     [Theta^0, Theta^1, ... Theta^6] * M * [Psi^0, Psi^1, Psi^2]' = 0

%%
% *Author:* Vladimir Grechka 1998 2012

%% 
% *Comment:*
%
% * Function |setMatVTI| is a modification of function |setMat| published in Grechka (2012)

%%
function [F, M] = setMatVTI(c11, c13, c33, c55)
%% Settings
[~, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

F = zeros(3,3);
M = zeros(7,3);

%% Construct the output matrices

% Assign nonzero elements of matrix F representing the Christoffel equation in the form
% [1, p1^2, p1^4] * F * [1, p3^2, p3^4]' = 0, where p1 and p3 are the horizontal and vertical
% components of the slowness vector 
F(1,1) = 1;
F(1,2) = -c33 - c55;
F(1,3) = c33*c55;
F(2,1) = -c11 - c55;
F(2,2) = c11*c33 - c13*(c13 + 2*c55);
F(3,1) = c11*c55;

if nargout == 1;   return;   end;

% Construct matrix M that relates Psi = tan(group angle) to Theta = tan(phase angle) as
% [Theta^0, Theta^1, ... Theta^6] * M * [Psi^0, Psi^1, Psi^2]' = 0
M(1,3) = F(1,3)*(-F(1,2)^2 + 4*F(1,1)*F(1,3));
M(2,2) = F(2,2)*( F(1,2)^2 - 4*F(1,1)*F(1,3));
M(3,1) = F(1,3)*F(2,1)^2 + F(2,2)*(-(F(1,2)*F(2,1)) + F(1,1)*F(2,2));
M(3,3) = -2*F(1,2)*F(1,3)*F(2,1) + 4*F(1,1)*F(1,3)*F(2,2);
M(4,2) = 2*(-(F(1,1)*F(2,2)^2) + F(1,2)^2*F(3,1) + F(1,3)*(F(2,1)^2 - 4*F(1,1)*F(3,1)));
M(5,1) = -2*F(1,2)*F(2,1)*F(3,1) + 4*F(1,1)*F(2,2)*F(3,1);
M(5,3) = -(F(1,2)*F(2,1)*F(2,2)) + F(1,1)*F(2,2)^2 + F(1,2)^2*F(3,1);
M(6,2) = F(2,2)*( F(2,1)^2 - 4*F(1,1)*F(3,1));
M(7,1) = F(3,1)*(-F(2,1)^2 + 4*F(1,1)*F(3,1));

end  % of the function