%% *timeLSQ*
% art3D-based LSQ function for locating events simultaneously with anisotropic velocity model 

%%
% *Input:*

% param     - vector of unknowns defined in function 'setUnknowns'
% printFlag - [scalar] flag, which should be equal to 1 to enable printing of results at 
%             every iteration

%%
% *Output:*

% F         - vector of traveltime differences
% J         - Frechet-derivative matrix (d F)/(d param) 

%%
% *Author:* Vladimir Grechka 2012 2013

%% 
% *Known issues:* 
%
% * check whether |reduceDerivMat| for the constraints is working correctly  

%%
function [F, J] = timeLSQ(param, printFlag)
%% Settings 
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

global units artFlags 
global dim data unknowns inversion modelCurrent artCurrent

%% Unpack the global structures into local variables

% The array sizes  
noPerf  = sum(modelCurrent.indexPerf);
noEvnt  = sum(modelCurrent.indexEvnt);
noSou   = noPerf + noEvnt;
noRec   = size(modelCurrent.xRec, 2);
noWave  = size(modelCurrent.waveCode, 2);
noTime  = noSou*noRec*noWave;
noParam = length(param);
noLayer = modelCurrent.noLayer;

%% Setup the output arrays
F = zeros((noTime + 10*noLayer), 1);
J = zeros((noTime + 10*noLayer), noParam);

% Set the flag for Jacobian calculation
artFlags.flagDT = 0;
if nargout > 1;  artFlags.flagDT = 1;  end;

%% Build current velocity model for ray tracing

% Stiffnesses
%noUnknCij = length(find(unknowns.CijInd == 0));
noUnknCij = numel(find(unknowns.CijInd == 0));
[modelCurrent.Cij.global, modelCurrent.Cij.symmetry, ...
    modelCurrent.rotation.angles, modelCurrent.rotation.matrix, ~, ~, ~, ~] = ...
    getCij(param(1:noUnknCij), modelCurrent.anisType, modelCurrent.anisParam, unknowns.CijInd, 0);
countUnkn = noUnknCij;

% Coordinates of the sources 
for isou = 1:noSou
    if unknowns.perfInd(isou) == 0
        % Unknown source location
        x0 = param((countUnkn + 1) : (countUnkn + dim));
        if dim == 2 
            modelCurrent.xSou(:,isou) = ...
                [modelCurrent.xWell(1) + x0(1)*cos(data.azimuth(isou)); ...
                 modelCurrent.xWell(2) + x0(1)*sin(data.azimuth(isou)); x0(2)];
        else
            modelCurrent.xSou(:,isou) = x0';
        end;
        countUnkn = countUnkn + dim;
    elseif unknowns.perfInd(isou) == 1
        % Known perforation-shot location -- do nothing        
%        modelCurrent.xSou(:,i) = xSou(:,i);
    else
        fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
        fprintf('>>> Parameter unknPrfInd(%g) = %g is invalid \n', [isou, unknPrfInd(isou)]);
        fprintf('>>> Its legitimate values are 0 or 1 \n \n');
          error('>>> STOP');
    end;
end;

% Interfaces
noUnknInt = length(find(unknowns.interfaceInd == 0));
[modelCurrent.interface, modelCurrent.noInt] = ...
    getInterface(unknowns.interfaceInd, param((countUnkn + 1) : (countUnkn + noUnknInt)));
countUnkn = countUnkn + noUnknInt;
           
% Faults
noUnknFlt = length(find(unknowns.faultInd == 0));
[modelCurrent.fault, modelCurrent.noFlt] = ...
    getFault(unknowns.faultInd, param((countUnkn + 1) : (countUnkn + noUnknFlt)));
countUnkn = countUnkn + noUnknFlt;

% Origin times
modelCurrent.tau = getTau(unknowns.tauInd, param((countUnkn + 1) : length(param)));

%% Compute traveltimes of the events and perforation shots
artCurrent = art3D(modelCurrent, artFlags, data.time);

%% Reduce the Frechet-derivative matrix to fit vector 'param'
if nargout > 1
    frechet = reduceDerivMat(artCurrent.dtdm, dim, modelCurrent, unknowns); 
end;
    
%% Weight the data  
weights = [1, 1];   % assign the weights to perf shots and microseismic events
% if noPerf ~= 0 && noEvnt ~= 0
%     weights = 0.5*(noPerf+noEvnt)*[1/noPerf, 1/noEvnt];  % equalize the relative weights of 
% elseif noPerf == 0                                       % all perfs and all events
%     weights = [0, 1];                  
% elseif noEvnt == 0;   
%     weights = [1, 0];   
% end;      

%% Build the objective function and its Jacobian 
count = 0;
for isou = 1:noSou
    % Assign the weights
    if modelCurrent.indexPerf(isou) == 1
        weight = weights(1);  % the weight for perf shots
    else
        weight = weights(2);  % the weight for microseismic events
    end;    

    for irec = 1:noRec
        for iwave = 1:noWave
            itime = noRec*noWave*(isou - 1) + noWave*(irec - 1) + iwave; 
            if isnan(artCurrent.time(itime,1)) == 1 || isnan(data.time(itime,1)) == 1
                % Either data point does not exist or time cannot be computed in the current model
                % --> keep F(itime,1) = 0 and J(itime,:) = 0;
                count = count + 1;
                data.time(itime) = NaN;  
                artCurrent.time(itime,1) = NaN;   artCurrent.traj(:,:,itime) = NaN;    
                artCurrent.n(:,:,itime)  = NaN;   artCurrent.r(:,:,itime) = NaN;   
                artCurrent.U(:,:,itime)  = NaN;   artCurrent.Kp(:,itime)  = NaN; 
                artCurrent.Vph(:,itime)  = NaN;   artCurrent.Vgr(:,itime) = NaN;          
            else
                artCurrent.time(itime,1) = artCurrent.time(itime) + modelCurrent.tau(isou);
                F(itime,1) = weight*(artCurrent.time(itime,1) - data.time(itime));
                if nargout > 1
                    J(itime,:) = weight*frechet(itime,:);
                end;                
            end;
        end;    % of loop over iwave
    end;    % of loop over irec
end;    % of loop over isou

inversion.timeMisfit = F;    inversion.dtdm = J;
inversion.timeRMS = norm(F)/sqrt(noTime - count);
inversion.percMissTraj = 100*count/noTime;

%% Penalize violatiion of the elastic stability conditions and the velocity constraints
factor = 1.e+5;  % enforce the constraints
if noUnknCij > 0
    [constrFun, constrJac] = ...
        setConstraints(modelCurrent.Cij.global, inversion.velBound, modelCurrent.wellVec, ...
        factor*inversion.timeRMS*[1,1]);
   
    % Reduce the size of the derivative matrix of constraints to the Cij portion of the parameter vector
    unknownsPenalty = unknowns;
    unknownsPenalty.interfaceInd = [];    unknownsPenalty.faultInd = [];
    unknownsPenalty.perfInd = [];         unknownsPenalty.tauInd = [];
    constrFrechet = reduceDerivMat(constrJac, dim, modelCurrent, unknownsPenalty); 
    F((noTime + 1) : (noTime + 10*noLayer), 1) = constrFun;
    J((noTime + 1) : (noTime + 10*noLayer), 1 : noUnknCij) = constrFrechet;
    inversion.penalty = norm(constrFun);
else
    inversion.penalty = 0;    % the constraints are satisfied: there is nothing to adjust   
end;

%% Apply the equality and inequality constraints
[Fcon, Jcon] = setParamConstraints(param);
F = cat(1, F, factor*inversion.timeRMS*Fcon);
J = cat(1, J, factor*inversion.timeRMS*Jcon);

%% Print the results at current iteration 
if printFlag == 1
    fprintf('>>> timeRMS (ms) = %13.6e,  penalty = %12.5e,  percMissTraj = %5.2f \n', ...
            u2u(inversion.timeRMS, units.time.processing, 'ms'), ...
            inversion.penalty, inversion.percMissTraj);
     disp('>>> Selected parameters:');
% %    disp(param([1 3 5 7  9])./param([2 4 6 8 10])); 
% %    disp(u2u(param([1:8, 14, 15, 21, 22, 28, 29]), 'm', 'kft'));  
     disp(param(1:10));  
%     disp(u2u(param([4, 5, 9, 10, 14, 15]), 'rad', 'deg'));  
end;

% display('>>> PAUSE');  pause;

end    % of the function