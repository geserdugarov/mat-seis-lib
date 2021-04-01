%%
% Generate traveltime data for differential event locations through NMO cylinders

clear variables;
close all;
warning off;  beep off;  format short;

% Determine the script name 
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));    
fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
fprintf('    is running... \n');

localMatlabFolder1 = '_VG Private M-Functions';  addpath(genpath(localMatlabFolder1));
localMatlabFolder2 = 'Local M-Functions';        addpath(genpath(localMatlabFolder2));

% setParPool;

%% Definitions
global modelflag sourecflag;
% see setCij
modelflag  = 7;
% 1 - ground profile, 2 - VSP, 3 - 2D survey, 4 - from Lykhin
% 5 - Srednee Priob'e, 6 - Verhnechonskoe, 7 - Turia 1
sourecflag = 7;

% Output file
outputFileName = '2HTIreflection';
disp(outputFileName);

% Figure parameters
fig = 0;    saveFigFlag = [];     
figPosition = figureGrid([2, 2], 50);

% Units
[units] = setUnits('m', 'm', 'm', 's', 's', 's', 'rad', 'rad', 'deg', 'Hz', 'Hz', 'Hz');
                        
% Flags for ray tracing
artFlags.flagWAA = 'SA';    %   'SA' is another option %    
artFlags.flagDT  = 1;    artFlags.flagDU  = 0;

%% Build anisotropic model, input acquisition geometry and ray codes
if modelflag==1
    model.anisType = ['VTI';'VTI';'VTI'];
elseif modelflag==2
    model.anisType = ['VTI';'VTI'];
elseif modelflag==3
    model.anisType = ['ISO';'ISO';'ISO'];
elseif modelflag==4
    model.anisType = ['ISO';'VTI';'ISO'];
elseif modelflag==5
    model.anisType = ['ISO';'ISO';'ISO';'ISO';'ISO';'ISO';'ISO';'VTI'];
elseif modelflag==6
    model.anisType = ['ISO';'ISO';'ISO';'ISO';'ISO';'ISO';'ISO';'VTI';'ISO';'ISO'];
elseif modelflag==7
    model.anisType = ['ISO';'ISO';'ISO';'VTI';'VTI'];
else
    model.anisType = ['VTI'];
end
model = setModel(model.anisType);

%% Ray tracing
disp('>>> --------------------------------------');
disp('>>> art3D is running...');
tic
artOut = art3D(model, artFlags, []);
toc

% checking polarization
noRays = size(artOut.rcode,3);
temp = artOut.rcode(1,:,1);
temp(isnan(temp)) = [];
% all rays have the same segment numbers
noRaySeg = numel(temp); % number of ray segments
clear temp;

Uall = artOut.U(:,1:noRaySeg,1:noRays);
for raynum = 1:noRays
    for RTnum = 2:noRaySeg
        % comparing with previous ray segment
        Uold = Uall(:,RTnum-1,raynum);
        Ucur = Uall(:,RTnum,raynum);        
        if abs(acos(dot(Ucur,Uold)) - pi/2) < 1.0e-6 && raynum >= 2
            % for SH-wave comparing with previous ray
            Uold = Uall(:,RTnum,raynum-1);
            Ucur = Uall(:,RTnum,raynum);
            if acos(dot(Ucur,Uold)) > pi/2
                Uall(:,RTnum,raynum) = -Uall(:,RTnum,raynum);
            end
        elseif acos(dot(Ucur,Uold)) > pi/2 + 1.0e-6
            Uall(:,RTnum,raynum) = -Uall(:,RTnum,raynum);
        end
    end
end
artOut.U(:,1:noRaySeg,1:noRays) = Uall;

%% Plot the model and rays
fig = fig + 1;
setFigure(fig, 3, figPosition(fig, :), {'Model'; outputFileName}, ...
    [ 'East (', units.length.display, ')'], ...
    ['North (', units.length.display, ')'], ...
    ['Depth (', units.length.display, ')'], [], []);
plotSouRec(fig, model.xSou, model.xRec, units.length.processing, units.length.display, ...
           model.indexPerf, 16, [1 0 0], model.indexEvnt, 7, [0 0 0], 10, [0 0 0], 0, 0);
newAxis = axis;        
plotPlane(fig, newAxis, [model.interface, model.fault], ...
    units.length.processing, units.length.display, [0,0,0], 0.2);
plotRayTraj(fig, artOut.traj, units.length.processing, units.length.display, model.noSou, model.noRec, model.noWave)
view(-8, 8);   % view(3);
%return

%% Dynamic ray tracing
disp('>>> --------------------------------------');
disp('\>>> rayDynam is running...');
tic
dynamOut = rayDynam(model, artOut);
toc


pos = find(isnan(artOut.traj(1,:,1)),1) - 1;
offset = vecnorm(reshape(artOut.traj(:,pos,:) - artOut.traj(:,1,:),3,size(artOut.traj,3)));
angles = squeeze(rad2deg(acos((artOut.n(3,4,:)))));
sinAng2 = sin(deg2rad(angles)).^2;
temp = [repmat(0,raynum,1) offset' angles sinAng2 (dynamOut.RTcoef(4,:))'];
% temp = [(0:5:360)' repmat(1500,raynum,1) angles' sinAng2' (dynamOut.RTcoef(4,:))'];

%% Calculating seismograms
% disp('>>> --------------------------------------');
% disp('>>> calculating seismograms...');
% pos = find(isnan(artOut.traj(1,:,1)),1) - 1; % there are some NaN values, find last position
% [time,seisZ] = calcTracesZ((artOut.time)', reshape(artOut.U(:,pos-1,:),3,size(artOut.U,3)), ...
%                                         dynamOut.RTcoef, ones(1,numel(dynamOut.L)), ...
%                                         'ricker', 1.5, 0.002);
% 
% offset = vecnorm(reshape(artOut.traj(:,pos,:) - artOut.traj(:,1,:),3,size(artOut.traj,3)));
% maxval = max(max(abs(seisZ)))/5;
% clear pos;

% pickspath = 'manualPick_X.dat';
% figX = figure;
% temp = seisX;
% set(figX,'units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
% subplot('Position',[0.05 0.06 0.90 0.88]); hold on;
% if maxval > 0
%     tempData = repmat(1:size(temp,2),size(temp,1),1) + temp./repmat(maxval,1,size(temp,2))/2;
% else
%     tempData = repmat(1:size(temp,2),size(temp,1),1) + temp./2;
% end
% plot(repmat(time,1,size(temp,2)),tempData,'r');
% xlim([0,max(time)]);
% ylim([0,41]);
% LinePlotExplorer_polyfit(pickspath);
% set(gca,'Ytick',1:size(temp,2)); set(gca,'YtickLabel',offset);
% grid on; xlabel('times, s'); ylabel('offset, m');
% title('X component')
% 
% pickspath = 'manualPick_Y.dat';
% figY = figure;
% temp = seisY;
% set(figY,'units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
% subplot('Position',[0.05 0.06 0.90 0.88]); hold on;
% tempData = repmat(1:size(temp,2),size(temp,1),1) + temp./repmat(maxval,1,size(temp,2))/2;
% plot(repmat(time,1,size(temp,2)),tempData,'g');
% xlim([0,max(time)]);
% ylim([0,41]);
% LinePlotExplorer_polyfit(pickspath);
% set(gca,'Ytick',1:size(temp,2)); set(gca,'YtickLabel',offset);
% grid on; xlabel('times, s'); ylabel('offset, m');
% title('Y component')

% pickspath = 'manualPick_Z.dat';
% figZ = figure;
% temp = seisZ;
% set(figZ,'units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
% subplot('Position',[0.05 0.06 0.90 0.88]); hold on;
% tempData = repmat(1:size(temp,2),size(temp,1),1) + temp./repmat(maxval,1,size(temp,2))/2;
% plot(repmat(time,1,size(temp,2)),tempData,'b');
% xlim([0,max(time)]);
% ylim([0,41]);
% LinePlotExplorer_polyfit(pickspath);
% set(gca,'Ytick',1:size(temp,2)); set(gca,'YtickLabel',offset);
% grid on; xlabel('times, s'); ylabel('offset, m');
% title('Z component')
% 
% offset = vecnorm(reshape(artOut.traj(:,15,:)-artOut.traj(:,1,:),[3 size(artOut.traj,3)]));
% waveNorm = reshape(artOut.n(:,7,:),[3 size(artOut.n,3)]);
% waveOut  = reshape(artOut.n(:,14,:),[3 size(artOut.n,3)]);
% upLayerRefr = prod(dynamOut.RTcoef([1:6,8:13],:),1);
% realRefl = dynamOut.RTcoef(7,:);

% save(strcat('az060_off50-50-2000_PP.mat'),'time','seisX','seisY','seisZ','offset','waveNorm');
% save(strcat('data\VCH-az000-off50-50-2000-PP.mat'),'time','seisZ','offset','waveNorm','waveOut','upLayerRefr','realRefl');

% % Load the results of kinematic ray tracing
% display('>>> --------------------------------------');
% fprintf('>>> art3D results have been loaded from file ''%s'' \n', ...
%     fullfile('Outputs', [outputFileName, '.mat']));
% load(fullfile('Outputs', [outputFileName, '.mat']));

% %% Change the model
% model.anisType = repmat('ISO', model.noLayer, 1);    
% model = setModel(model.anisType);

%% Compute NMO-velocity cylinders
% display('>>> --------------------------------------');
% display('>>> artNMO is running...');
% nmo = artNMO(model, artOut);

%% Save the ray-tracing results
% save(fullfile('Outputs', [outputFileName, '.mat']), '-mat', ...
%     'model', 'units', 'artFlags', 'artOut', 'nmo'); 

%%
% rmpath(genpath(localMatlabFolder1));  % remove local folders from the path
% rmpath(genpath(localMatlabFolder2));  
fprintf('>>> Return from script ''%s'' \n \n', [thisFileName, '.m']);  return          
