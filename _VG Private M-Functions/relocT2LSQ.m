%% *relocT2LSQ*
% Objective function for the least-squares optimization of relative event locations
% based on the paraxial ray series for T^2

%%
% *Input:*

% param   - [4, 1] vector containing three spatail coordinates of the hypocenter of a slave event 
%           relative to a master event, and the origin time of the slave
% tMaster - [:, 1] array of times of the master event
% pMaster - [3, :] array of the slowness vectors evaluated at the master source location of the  
% UMaster - [3, 3, :] array of the NMO-velocity surfaces evaluated at the master source location 
% tSlave  - [:, 1] array of times of the slave event

%%  
% *Output:*

% F       - [:, 1] objective function to be minimized
% J       - [:, 4] Jacobian of the objective function  

%%
% *Author:* Vladimir Grechka 2014

%%
function [F, J] = relocT2LSQ(param, tMaster, pMaster, UMaster, tSlave)
%% Settings 
[~, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

global FF Jac numCondPCH

dxSou = param(1:3,1);  dtau = param(4,1);
F = NaN(size(tMaster));
J = NaN(size(tMaster,1), size(param,1));

%% Second-order traveltime extrapolation
for idata = 1:length(tMaster)
    tmp = kron(pMaster(:,idata)', pMaster(:,idata)) + UMaster(:,:,idata);
    t2Predicted = tMaster(idata, 1)^2 + 2*tMaster(idata, 1)*dot(pMaster(:,idata), dxSou) ...
                + dxSou'*tmp*dxSou; 
    F(idata,1) = t2Predicted - (tSlave(idata, 1) - dtau)^2;
    J(idata,:) = 2*[tMaster(idata, 1)*pMaster(:,idata)' + dxSou'*tmp, tSlave(idata, 1)];
end;

%% Outputs for further analysis (can be commented out)
% Since F is the difference of the SQUARED times, F = tm^2 - ts^2 ~ tm^2 - (tm^2 + 2*tm*dt)
FF = 0.5*norm(F./tMaster)/sqrt(length(tMaster));
Jac = J;

numCondPCH = NaN;
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
