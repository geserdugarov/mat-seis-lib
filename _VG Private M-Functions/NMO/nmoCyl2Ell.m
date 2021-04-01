%% *nmoCyl2Ell*
% Compute an NMO ellipse as a cross-section of an NMO cylinder by a plane

%%
% *Input:*
% Ucyl - [3, 3] matrix expressing the NMO cylinder
% nrm  - [3, 1] vector normal to the cross-section plane

%%
% *Output:*
% W    - [2, 2] symmetric matrix representing the NMO ellipse

%%
% Author: Vladimir Grechka 1998, 2014

%%
% *Comments:*
%
% * Original Matlab version, called 'cyl2ell', is a part of the 'ART' package freely  
%   distributed by the CWP, CSM

%%
function [W] = nmoCyl2Ell(Ucyl, nrm)
%% Settings  
[~, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

W = zeros(2,2);

%% Construct vectors normal to 'nrm' that lie in the plane of the sought NMO ellipse
[azm, pol90, ~] = cart2sph(nrm(1), nrm(2), nrm(3));
pol = pi/2 - pol90;
b(:,1) = [cos(pol)*cos(azm), cos(pol)*sin(azm), -sin(pol)]';
b(:,2) = [-sin(azm), cos(azm), 0]';

%b
%[b, ~, ~] = dVecdAngle(nrm)

%% Build a cross-section of the NMO cylinder (equations D7 and D8 in Grechka and Tsvankin, 2002)
for i=1:2
    for j=1:2
        for k=1:3
            for m=1:3
                W(i,j) = W(i,j) + (b(k,i)*b(m,j) + b(k,j)*b(m,i))*Ucyl(k,m)/2;
            end;
        end;
    end;
end;

end    % of the function

