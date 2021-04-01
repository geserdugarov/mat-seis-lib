%% *rec2geoLSQ3D*
% Objective function for the least-squares optimization of receiver orientations in 3D

%%
% *Input:*

% RotVec  - [9, 1] vector representation of the [3, 3] rotation matrix

%%  
% *Output:*

% F     - [scalar] objective function to be minimized

%%
% *Author:* Vladimir Grechka 2014

%%
function [F, J] = timeDiff(param, tMaster, pMaster, UMaster, tSlave)
%% Settings 
[~, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

xSlave = param(1:2,1);  
tauSlave = param(3,1);
F = NaN(size(tMaster));
J = NaN(size(tMaster,1), size(param,1));

for idata = 1:length(tMaster)
    F(idata,1) = ...
        tauSlave + tMaster(idata, 1) ...
                 + dot(pMaster(:,idata), xSlave) ...
                 + xSlave'*UMaster(:,:,idata)*xSlave/(2*tMaster(idata, 1)) ...
                 - tSlave(idata, 1); 
%    J(idata,:) = [pMaster(idata,:), 1];
    J(idata,:) = [pMaster(:,idata)' + xSlave'*UMaster(:,:,idata)/tMaster(idata, 1), 1];
end;
        

end  % of the script
