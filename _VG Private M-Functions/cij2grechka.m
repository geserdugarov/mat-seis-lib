%% *cij2grechka*
% Calculate Grechka's (2000) anisotropy parameters and anellipticity coefficients for a given 
% (presumably) monoclinic stiffness matrix 

%% 
% *Input:*  

% Cij - [6, 6] stiffness matrix 

%%
% *Output:*

% Ani - [1, 15] vector of anisotropy parameters 
%       [Vp0, Vs0, epsilon1, epsilon2, delta1, delta2, delta3, gamma1, gamma2, 
%        zeta1, zeta2, zeta3, eta1, eta2, eta3]

%%
% *Author:* Vladimir Grechka 2000 2012

%%
function [Ani] = cij2grechka(Cij)
%% Settings and defaults
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

if isequal(size(Cij), [6, 6]) == 0
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Incorrect size(Cij) = [%g, %g] \n \n', size(Cij));   
      error('>>> STOP');
end;

%% Definitions of the anisotropy parameters
Vp0 = sqrt(Cij(3,3));   
Vs0 = sqrt(Cij(5,5));
epsilon2 = (Cij(1,1) - Cij(3,3))/(2*Cij(3,3));    
epsilon1 = (Cij(2,2) - Cij(3,3))/(2*Cij(3,3));
delta2 = ( (Cij(1,3) + Cij(5,5))^2 - (Cij(3,3) - Cij(5,5))^2 )/(2*Cij(3,3)*(Cij(3,3) - Cij(5,5)));
delta1 = ( (Cij(2,3) + Cij(4,4))^2 - (Cij(3,3) - Cij(4,4))^2 )/(2*Cij(3,3)*(Cij(3,3) - Cij(4,4)));
delta3 = ( (Cij(1,2) + Cij(6,6))^2 - (Cij(1,1) - Cij(6,6))^2 )/(2*Cij(1,1)*(Cij(1,1) - Cij(6,6)));
gamma2 = (Cij(6,6) - Cij(4,4))/(2*Cij(4,4));
gamma1 = (Cij(6,6) - Cij(5,5))/(2*Cij(5,5));
zeta1  = (Cij(1,6) - Cij(3,6))/(2*Cij(3,3));
zeta2  = (Cij(2,6) - Cij(3,6))/(2*Cij(3,3));
zeta3  = Cij(3,6)/Cij(3,3);
eta1   = (epsilon1 - delta1)/(1 + 2*delta1);
eta2   = (epsilon2 - delta2)/(1 + 2*delta2);
eta3   = (epsilon1 - epsilon2 - delta3*(1 + 2*epsilon2))/((1 + 2*epsilon2)*(1 + 2*delta3));

Ani = [Vp0, Vs0, epsilon1, epsilon2, delta1, delta2, delta3, gamma1, gamma2, ...
       zeta1, zeta2, zeta3, eta1, eta2, eta3];

end    % of the function