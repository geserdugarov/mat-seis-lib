%% *isSnell*
% Check whether Snell's law at model interfaces is satisfied with the prescribed tolerance

%%
% *Input:*

% p       - [3, :] array of slowness vectors along the ray trajectory
% nrmInt  - [3, :] array of normals to the interfaces intersected by the ray
% tol     - [scalar] tolerance

%%
% *Output:*

% flagInt - [size(nrmInt, 2), 1] array whose elements are equal to 1 if Snell's law is satisfied
%           at an interface with tolerance 'tol' and 0 otherwise
% flagCum - [scalar] equal to 1 if Snell's law is satisfied for all interfaces and 0 otherwise

%%
% *Author:* Vladimir Grechka 2014

%%
function [flagInt, flagCum] = isSnell(p, nrmInt, tol)
%% Settings and defaults
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
% narginTrue = nargin(thisFileName);   
% narginchk(narginTrue, narginTrue);

if nargin < 3;  tol = 1.e-4;  end;   % default for the tolerance
flagInt = ones(size(nrmInt, 2), 1);
flagCum = 1;   

%% Check Snell's law
if size(p, 2) == 1
    flagInt = [];
    return;
    
else
    for iseg = 1:size(nrmInt, 2)
        tmp1 = cross(p(:, iseg),     nrmInt(:, iseg)) - ...
               cross(p(:, iseg + 1), nrmInt(:, iseg));
        tmp2 = (norm(p(:, iseg)) + norm(p(:, iseg + 1)))/2;
        if norm(tmp1) < tol*norm(tmp2)
            flagInt(iseg, 1) = 1;
        else
            flagInt(iseg, 1) = 0;
            flagCum = 0;
        end;
%        fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
%        fprintf('>>> iseg = %g,  diff = %g, aver = %g,  diff/aver = %g \n',  ...
%            iseg, norm(tmp1), norm(tmp2), norm(tmp1)/norm(tmp2) );
    end;
end;
        
end    % of the function

%%