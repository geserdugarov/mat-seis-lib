%% *vec2axis*
% Construct rotation matrix |R| that maps |R*vec1| onto axis |x(j)| and places |R*vec2| in the 
% plane orthogonal to |x(k)|

%%
% *Input:*

% vec1 - [3, 1] first input vector that becomes the coordinate axis x(j) after rotation
% j    - [scalar] equal to 1, 2, or 3 and indicating the number of coordinate axis pointing 
%        along R*vec1
% vec2 - [3, 1] second input vector placed in the plane orthogonal to x(k) after rotation
% k    - [scalar] equal to 1, 2, or 3 (k ~= j) and indicating the number of coordinate axis  
%        orthogonal to R*vec2; if x = R*vec2, than x(k) = 0

%%
% *Output:*

% R    - [3, 3] unitary matrix

%%
% *Author:* Vladimir Grechka 1998 2012 

%%
function [R] = vec2axis(vec1, j, vec2, k)
%% Settings and checks
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

if isempty(vec1) == 1
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Variable ''vec1'' is undefined \n');
      error('>>> STOP');
end;

if isempty(vec2) == 1;    vec2 = vec1;    end;

if isempty(j) == 1
    j = 1;
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Coordinate axis ''j'' is undefined -- defaulted to %g \n', j);
end;

if isempty(k) == 1
    tmp = [2, 3, 1];    k = tmp(j);
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Coordinate axis ''k'' is undefined -- defaulted to %g \n', k);
end;

if j == k
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Improper equality of the coordinate axes: j = %g and k = %g \n', [j, k]);
      error('>>> STOP');
end;

tol  = 1.e-8;   
vec1 = vec1/norm(vec1);     % normalize the input vectors
vec2 = vec2/norm(vec2);

%% Construct the rotation matrix
x1 = vec1;   crosspr = cross(vec2, x1);   
if norm(crosspr) < tol
    dv = dVecdAngle(vec1);  % avoid division by zero 
    x2 = dv(:,2)/norm(dv(:,2));
else
    x2 = crosspr/norm(crosspr);
end;
x3 = cross(x1, x2);
R = [x1, x2, x3]';          % this matrix R corresponds to j = 1, k = 2

%% Switch the coordinate directions  
switch j
    case 1;
        switch k
            case 2;         % do nothing                % j = 1, k = 2
            case 3;   R = [R(1,:); R(3,:); R(2,:)];     % j = 1, k = 3
        end;
    case 2;    
        switch k
            case 1;   R = [R(2,:); R(1,:); R(3,:)];     % j = 2, k = 1
            case 3;   R = [R(3,:); R(1,:); R(2,:)];     % j = 2, k = 3
        end;
    case 3;    
        switch k
            case 1;   R = [R(2,:); R(3,:); R(1,:)];     % j = 3, k = 1
            case 2;   R = [R(3,:); R(2,:); R(1,:)];     % j = 3, k = 2
        end;
end;

end  % of the function
