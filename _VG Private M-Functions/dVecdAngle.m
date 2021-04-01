%% *dVecdAngle*
% Calculate the first- and second-order derivatives of the unit vector with respect to its
% spherical angles  

%%
% *Input:*

% vec - [3, 1] unit vector given by
%       vec = [sin(beta(1))*cos(beta(2)), sin(beta(1))*sin(beta(2)), cos(beta(1))]' 

%%
% *Output:*

%  d1  - [3, 2] array of derivatives d vec/d beta(i)
%        d1(:,i) = d vec/d beta(i) 
%  spn = sqrt(vec(1)^2 + vec(2)^2) 
%  d2  - [3, 3] array of second-order derivatives 
%        d2(:,i+j-1) = d^2 vec/(d beta(i) d beta(j))

%%
% *Author:* Vladimir Grechka 1989 1998 2012
%
% * Fortran version is published in Obolentseva and Grechka (1989)

%%
function [d1, spn, d2] = dVecdAngle(vec)
%% Settings  
[~, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

%% Normalize variable 'vec'
vec = vec/norm(vec);                                                
spn = sqrt(vec(1)^2 + vec(2)^2);

%% Avoid the coordinate singularity at beta(1) = 0
tol1 = 1.e-8;  tol2 = 1.e-6;                                        % set the tolerances
if spn < tol1
    vec(1) = tol2;
    vec = vec/norm(vec);                                            % renormalize variable vec
    spn = sqrt(vec(1)^2 + vec(2)^2);                                % recompute variable spn
end;

%% The first-order derivatives
d1(:,1) = [ vec(1)*vec(3)/spn,  vec(2)*vec(3)/spn, -spn]';          
d1(:,2) = [-vec(2),             vec(1),               0]';

%% The second-order derivatives
if nargout > 2
    d2(:,1) =  -vec;                                                
    d2(:,2) = [-vec(2)*vec(3)/spn,  vec(1)*vec(3)/spn, 0]';
    d2(:,3) = [-vec(1),            -vec(2),            0]';
end;

end    % of the function