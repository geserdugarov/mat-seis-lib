%% *relocLSQ*
% Objective function for the least-squares optimization of relative event locations

%%
% *Input:*

% param    - [4, 1] vector containing three spatail coordinates of the hypocenter of a slave event 
%            relative to a master event, and the origin time of the slave
% tMaster  - [:, m] array of times of the m-th master events
% pMaster  - [3, :, m] array of the slowness vectors evaluated at the master sources location of the  
% UMaster  - [3, 3, :, m] array of the NMO-velocity surfaces evaluated at the master source locations
% tSlave   - [:, 1] array of times of the slave event
% dxMaster - [3, m] array of differences between the first and m-th master event coordinates

%%  
% *Output:*

% F       - [:, 1] objective function to be minimized
% J       - [:, 4] Jacobian of the objective function  

%%
% *Author:* Vladimir Grechka 2014

%%
function [F, J] = relocLSQ(param, tMaster, pMaster, UMaster, tSlave, dxMaster)
%% Settings 
[~, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);
narginchk(narginTrue, narginTrue);

global FF Jac numCondPCH

dxSou = param(1:3,1);  dtau = param(4,1);
F = NaN(size(tMaster,1)*size(tMaster,2), 1);
J = NaN(size(tMaster,1)*size(tMaster,2), size(param,1));

%% Second-order traveltime extrapolation
for mdata = 1:size(tMaster,1)
    for idata = 1:size(tMaster,2)
    tPredicted = tMaster(mdata,idata) + dot(pMaster(:,idata,mdata), dxSou) ...
               + (1/2)*dxSou'*UMaster(:,:,idata,mdata)*dxSou/(tMaster(mdata,idata, 1));
	fewMasters = (dxMaster(:,mdata))'*UMaster(:,:,idata,mdata)*dxSou/(tMaster(mdata,idata, 1)) ...
               + dot(pMaster(:,idata,mdata), dxMaster(:,mdata)) ...
               + (1/2)*(dxMaster(:,mdata))'*UMaster(:,:,idata,mdata)*dxMaster(:,mdata)/(tMaster(mdata,idata, 1));
    F((mdata-1)*size(tMaster,2)+idata,1) = tPredicted + fewMasters - (tSlave(idata, 1) - dtau);
    J((mdata-1)*size(tMaster,2)+idata,:) = [pMaster(:,idata,mdata)' + dxSou'*UMaster(:,:,idata,mdata)/tMaster(mdata,idata), 1];
    end;
end;

%% Outputs for further analysis (can be commented out)
% FF = norm(F)/sqrt(length(tMaster));
% Jac = J;
% 
% Mat = Jac;  tolMax = 1.e+10;  % tolerance intended to remove the amplification of round-off errors
% scaleOut = ones(1,size(Mat,2));
% for i=1:size(Mat,2)
%     scaleOut(1,i) = 1/norm(Mat(:,i));
%     if abs(scaleOut(1,i)) > tolMax
%         scaleOut(1,i) = tolMax;
%     end;
% end;
% MatOut = Mat*diag(scaleOut);  numCondPCH = cond(MatOut);
    
end  % of the script
