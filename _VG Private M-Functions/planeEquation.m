%% *planeEquation*
% Calculate coefficients describing planar interface in the form |z = A(1)*x + A(2)*y + A(3)|  

%%
% *Input:*

% plane - [5, :] array specified by functions 'setInterface' or 'setFault'

%%
% *Output:*

% A     - [3, :] array of z-coefficients in equation z = A(1)*x + A(2)*y + A(3) 
% N     - [3, :] array of the unit plane normals 

%%
% *Author:* Vladimir Grechka 2012

%%
function [A, N] = planeEquation(plane)
%% Settings and checks
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

noPln = size(plane, 2);
A = zeros(3, noPln);
N = zeros(3, noPln);
if isempty(plane) == 1;    return;   end;

%% Construct the equation describing a plane
for i = 1:noPln
    nrmZ  = 1 - (plane(4,i)^2 + plane(5,i)^2);
    if nrmZ <= 0 
        fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
        fprintf('>>> Vertical component of the plane normal is equal to \n');
        display(sqrt(nrmZ));
        fprintf('>>> Such a plane cannot exist \n \n');
          error('>>> STOP');
    else
        N(  :,i) = [plane(4,i); plane(5,i); sqrt(nrmZ)];
        A(  3,i) = dot(plane(1:3,i), N(:,i))/N(3,i);
        A(1:2,i) = -N(1:2,i)/N(3,i);
    end;
end;

end    % of the function