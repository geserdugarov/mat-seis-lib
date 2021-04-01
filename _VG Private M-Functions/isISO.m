%% *isISO*
% Check whether input stiffness matrix Cij is isotropic within a given tolerance

%%
% *Input:*

% Cij  - [6, 6] stiffness matrix
% tol  - [scalar] tolerance of proximity to isotropy

%%
% *Output:*

% flag - [scalar] equal to 1 if Cij is isotropic within given tolerance and to 0 otherwise

%%
% *Author:* Vladimir Grechka 1998 2012 - 2014
%
% * Fortran version is published in Obolentseva and Grechka (1989)

%%
function [flag] = isISO(Cij, tol)
%% Settings and defaults

% [~, thisFileName, ~] = fileparts(mfilename('fullpath'));
% narginTrue = nargin(thisFileName);   
% narginchk(narginTrue, narginTrue);

if nargin < 2  ||  isempty(tol) == 1
    tol = 1.e-6;                        % default for 'tol'
end        

%% The best-fit isotropic approximation of Cij (Fedorov, 1968) 
sum0 = Cij(1,1) + Cij(2,2) + Cij(3,3);
sum1 = sum0 + 2*(Cij(1,2) + Cij(1,3) + Cij(2,3));
sum2 = sum0 + 2*(Cij(4,4) + Cij(5,5) + Cij(6,6));
lambda = (2*sum1 - sum2)/15;   mu = (3*sum2 - sum1)/30;  
Ciso = thomsen2cij([sqrt(lambda + 2*mu), sqrt(mu), 0, 0, 0]); 

dC = norm(Cij - Ciso)/norm(Cij);        % relative difference between Cij and Ciso

if dC < tol
    flag = 1;                           % isotropy
else
    flag = 0;                           % anisotropy
end

end    % of the function