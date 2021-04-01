%% *isStable*
% Check whether given stiffness matrix is real-valued, symmetric, and stable  

%%
% *Input:*

% Cij    - [6, 6] stiffness matrix
% msg    - [scalar] turning on (msg = 1, default) and off (msg = 0) a verbose warning message

%%
% *Output:*

% flag   - [scalar] equal to 1 if Cij complies with all conditions and to 0 otherwise
% detCij - [6, 1] array of six determinants of the principal minors of matrix Cij

%%
% *Author:* Vladimir Grechka 1998 2012 - 2014

%%
function [flag, detCij] = isStable(Cij, msg)
%% Settings and defaults
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
% narginTrue = nargin(thisFileName);   
% narginchk(narginTrue, narginTrue);

if nargin < 2  ||  msg ~= 0;   msg = 1;   end;   % default for printing verbose messages
flag = 1;   tol = 1.e-12;

%% Calculate the minors of Cij
detCij = ones(6, 1);    sgnCij = ones(6, 1);
for i=1:6
    detCij(i,1) = det(Cij(1:i,1:i));
    sgnCij(i,1) = sign(detCij(i,1)); 
end;

%% Set the flag
if isreal(Cij) == 0  ||  max(max(abs((Cij + Cij')/2 - Cij)))/norm(Cij) > tol  ||  sum(sgnCij) < 6
    flag = 0;
end;

%% Check the Poisson's ratio for isotropy
if isISO(Cij, tol) == 1
    nu = Cij(1,3)/(2*(Cij(1,3) + Cij(4,4)));
    if nu < 0  &&  msg == 1
        fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
        fprintf('>>> WARNING: Negative Poisson''s ratio = %g \n', nu); 
        fprintf('             for Cij = \n');
        disp(Cij);
    end;
end;

%% Issue the warning message
if flag == 0  &&  msg == 1
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Stiffness tensor Cij is either asymmetric, or complex-valued, \n');
    fprintf('    or violates the elastic stability conditions \n');
    display(Cij)
    fprintf('\n>>> det = [%g, %g, %g, %g, %g, %g] \n', detCij);
    fprintf('\n>>> PAUSE \n \n');   pause
end;

end    % of the function