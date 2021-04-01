%% *planeDepth*
% Calculate the depth of a planar interface  

%%
% *Input:*

% plane - [5, :] array describing planar interface or fault 
%         (see functions 'setInterface' and 'setFault')
% xy    - [2, 1] vector of the lateral coordinates at which the depths are to be calculated 
% index - [1, :] array of the plane numbers whose depths are to be calculated

%%
% *Output:*

% depth - [1 : length(index)] array of the calculated depths

%%
% *Author:* Vladimir Grechka 2012

%%
function [depth] = planeDepth(plane, xy, index)
%% Settings 
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

noPln = size(plane,2);
depth = NaN(noPln,1);

%% Checks
if isempty(plane) == 1;    depth = NaN(size(index));    return;   end;

% Consistency check
if min(index) <= 0  ||  max(index) > noPln
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    display(index);
    fprintf('>>> Erroneous array of plane numbers \n \n');
      error('>>> STOP');
end;

%% Calculate the plane depths  
for i = 1:length(index)
    xyPln = plane(1:2,index(i));  
    zPln  = plane(3,index(i));
    nrmXY = plane(4:5,index(i));
    nrmZ  = 1 - (plane(4, i)^2 + plane(5, i)^2);
    if nrmZ <= 0 
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
        fprintf('>>> Vertical component of the plane normal is equal to \n');
        sqrt(nrmZ)
        depth(index(i)) = NaN;
    else
        nrmZ = sqrt(nrmZ);
        depth(index(i)) = zPln - dot(nrmXY, (xy - xyPln))/nrmZ;
    end;
    
end;

end    % of the function