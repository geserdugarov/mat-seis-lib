clear variables;
close all;
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
fprintf('    is running... \n');

localMatlabFolder1 = 'FuncBase';  addpath(genpath(localMatlabFolder1));
localMatlabFolder1 = '_VG Private M-Functions';  addpath(genpath(localMatlabFolder1));

%% Открытие текстового файла для чтения данных;
fileID = fopen('forInverseData_PP_PS.txt','r');
if fileID == -1 
    error('File is not opened'); 
end 
temp = fscanf(fileID,'%f %f %f\n', [3 Inf]);
temp = temp';

angAzim  = temp(:,1);
angIncid = temp(:,2);
Epp      = temp(:,3);

waveNum = 1;

%% inverse problem
% synthetic
densR = 2200;   VpR = 4660;   VsR = 2330;
epsR = 0;   gamR = 0;   delR = 0;   axAzimR = 0;   axTiltR = 0;
densT = 2600;   axAzimT = 60;   axTiltT = 90;
VpV = 5920; VsV = 2900;
epsV = -0.039;  deltaV = -0.194; gamT = 0.253;
epsT = -epsV/(1+2*epsV);
VpT = VpV/sqrt(1+2*epsT);
VsT = VsV/sqrt(1+2*gamT);
f = 1-(VsT/VpT)^2;
delT = deltaV*(1+2*epsT)*(1+2*epsT/f)+2*epsT*(1+epsT/f);

% Turia
% densR = 2400;   VpR = 3400;   VsR = VpR*0.5;
% epsR = 0;   gamR = 0;   delR = 0;   axAzimR = 0;   axTiltR = 0;
% densT = 2500;   VpTv = 4300;   VsTv = VpTv*0.55;
% epsTv = -0.01;   gamTv = 0.01;   delTv = -0.01;   axAzimT = 0;   axTiltT = 90;
% epsT = -epsTv/(1+2*epsTv);  gamT = gamTv; %gamT = -gamTv/(1+2*gamTv);
% VpT = VpTv/sqrt(1+2*epsT); VsT = VsTv/sqrt(1+2*gamT);
% f = 1-VsT^2/VpT^2;
% delT = delTv*(1+2*epsT)*(1+2*epsT/f) + 2*epsT*(1+epsT/f);

CijR0 = thomsen2cij([VpR, VsR, epsR, delR, gamR]);
CijT0 = thomsen2cij([VpT, VsT, epsT, delT, gamT]);

rotAngleR = u2u([axAzimR, axTiltR, 0]', 'deg', 'rad');
rotAngleT = u2u([axAzimT, axTiltT, 0]', 'deg', 'rad');
% Apply rotation  
if ~isempty(find(rotAngleR ~= 0, 1))
    % Rotate the stiffness matrix
    rotMatrix = loc2glb(rotAngleR);
    CijR = bond(CijR0, rotMatrix);
else
    CijR = CijR0;
end
if ~isempty(find(rotAngleT ~= 0, 1)) 
    % Rotate the stiffness matrix
    rotMatrix = loc2glb(rotAngleT);
    CijT = bond(CijT0, rotMatrix);
else
    CijT = CijT0;
end


% RTout = (RTcoefModifForInverse(CijR,CijT,densR,densT,deg2rad((0:1:360)'),deg2rad(repmat(30,numel(0:1:360),1)),1,1))';
% RTout = (RTcoefModifForInverse(CijR,CijT,densR,densT,deg2rad(30),deg2rad(30),1,1))';
% plot(RTout)

RTout = (RTcoefModifForInverse(CijR,CijT,densR,densT,deg2rad(repmat(0,numel(0:1:60),1)),deg2rad((0:1:90)'),1,1))';
plot(real(RTout))

tic
RTout = (RTcoefModifForInverse(CijR,CijT,densR,densT,deg2rad(angAzim),deg2rad(angIncid),1,1))';
toc

fig = figure;
set(fig,'units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
subplot('Position',[0.05 0.09 0.90 0.88]); hold on; grid on;
multi = RTout(1,1)/sqrt(Epp(1));
plot(abs(RTout),'*c')
plot(sqrt(Epp)*multi,'xk');

return


densR = 2200;   VpR = 4660;   VsR = 2330;
epsR = 0;   gamR = 0;   delR = 0;   axAzimR = 0;   axTiltR = 0;

densT = 2600;   axAzimT = 60;   axTiltT = 90;

xReal = [5683; 2367.1; 0.0423; -0.1445; 0.253012];
CijR0 = thomsen2cij([VpR, VsR, epsR, gamR, delR]);

CijT0 = @(x) thomsen2cij([x(1), x(2), x(3), x(4), x(5)]);
x0 = [6000; 3000; 0.01; 0.01; 0.01];
% CijT0 = @(x) thomsen2cij([x(1), x(2), xReal(3), xReal(4), xReal(5)]);
% x0 = [6000; 3000];

rotAngleR = u2u([axAzimR, axTiltR, 0]', 'deg', 'rad');
rotAngleT = u2u([axAzimT, axTiltT, 0]', 'deg', 'rad');
% Apply rotation  
if ~isempty(find(rotAngleR ~= 0, 1))
    % Rotate the stiffness matrix
    rotMatrix = loc2glb(rotAngleR);
    CijR = bond(CijR0, rotMatrix);
else
    CijR = CijR0;
end
if ~isempty(find(rotAngleT ~= 0, 1)) 
    % Rotate the stiffness matrix
    rotMatrix = loc2glb(rotAngleT);
    CijT = @(x) bond(CijT0(x), rotMatrix);
else
    CijT = @(x) CijT0(x);
end

tic
% funcMin = @(x) sqrt(sum( (RTcoefModifForInverse_IsoHTI(CijR,CijT(x),densR,densT,waveNormInverse(:,:,1),[0;0;1],1,1) - sqrt(enerSignalNoiseInverse(:,1)')*multi).^2 )) + ...
%                sqrt(sum( (RTcoefModifForInverse_IsoHTI(CijR,CijT(x),densR,densT,waveNormInverse(:,:,2),[0;0;1],1,2) - sqrt(enerSignalNoiseInverse(:,2)')*multi).^2 ));
funcMin = @(x) sqrt(sum( (RTcoefModifForInverse(CijR,CijT(x),densR,densT,waveNormInverse(:,:,1),[0;0;1],1,1) - sqrt(enerSignalNoiseInverse(:,1)')*multi).^2 )) + ...
               sqrt(sum( (RTcoefModifForInverse(CijR,CijT(x),densR,densT,waveNormInverse(:,:,2),[0;0;1],1,2) - sqrt(enerSignalNoiseInverse(:,2)')*multi).^2 ));
options = optimset('OutputFcn',@outFunWithXval,'Display','iter');
[x,fval,exitflag,output] = fminsearch(funcMin,x0,options);
toc

return




densR = 2200;   VpR = 4660;   VsR = 2330;
epsR = 0;   gamR = 0;   delR = 0;   axAzimR = 0;   axTiltR = 0;

densT = 2600;   axAzimT = 60;   axTiltT = 90;

xReal = [5683; 2367.1; 0.0423; -0.1445; 0.253012];
CijR0 = thomsen2cij([VpR, VsR, epsR, gamR, delR]);
CijT0 = @(x) thomsen2cij([x(1), x(2), xReal(3), xReal(4), xReal(5)]);

rotAngleR = u2u([axAzimR, axTiltR, 0]', 'deg', 'rad');
rotAngleT = u2u([axAzimT, axTiltT, 0]', 'deg', 'rad');
% Apply rotation  
if ~isempty(find(rotAngleR ~= 0, 1))
    % Rotate the stiffness matrix
    rotMatrix = loc2glb(rotAngleR);
    CijR = bond(CijR0, rotMatrix);
else
    CijR = CijR0;
end
if ~isempty(find(rotAngleT ~= 0, 1)) 
    % Rotate the stiffness matrix
    rotMatrix = loc2glb(rotAngleT);
    CijT = @(x) bond(CijT0(x), rotMatrix);
else
    CijT = @(x) CijT0(x);
end

tic
funcMin = @(x) sqrt(sum((RTcoefModifForInverse(CijR,CijT(x),densR,densT,waveNormInverse(:,:,1),[0;0;1],1,1).^2*multi - enerSignalNoiseInverse(:,1)').^2)) + ...
               sqrt(sum((RTcoefModifForInverse(CijR,CijT(x),densR,densT,waveNormInverse(:,:,2),[0;0;1],1,2).^2*multi - enerSignalNoiseInverse(:,2)').^2));
valsVp = 4500:500:6500;
valsVs = 1500:500:3500;
xReal = [5683; 2367.1; 0.0423; -0.1445; 0.253012];
funcMinVals = zeros(numel(valsVp)*numel(valsVs),3);
count = 1;
for idVp = 1:numel(valsVp)
    for idVs = 1:numel(valsVs)
        fprintf(['\nStep ',num2str(count),' from ',num2str(size(funcNimVals,1))]);
        Vp = valsVp(idVp);
        Vs = valsVs(idVs);
        funcMinVals(count,:) = [Vp Vs funcMin([Vp; Vs])];
        count = count + 1;
    end
end
toc

X = reshape(funcMinVals(:,1),numel(valsVp),numel(valsVs));
Y = reshape(funcMinVals(:,2),numel(valsVp),numel(valsVs));
Z = sqrt(reshape(funcMinVals(:,3),numel(valsVp),numel(valsVs)));
figure;
surf(X,Y,Z);


fprintf('\n>>> Return from script ''%s'' \n \n', [thisFileName, '.m']);
clear thisFileName thisFolderName;
return


x = [5657.598; 2400.522; 0.076; -0.12; 0.278];
xReal = [5683; 2367.1; 0.0423; -0.1445; 0.253012];
(xReal - x)./xReal*100


plot(multi*RTcoefModifForInverse(CijR,CijT(x),densR,densT,waveNormInverse(:,:,1),[0;0;1],1,1).^2,'or')
plot(multi*RTcoefModifForInverse(CijR,CijT(x),densR,densT,waveNormInverse(:,:,2),[0;0;1],1,2).^2,'or')

fprintf('>>> Return from script ''%s'' \n \n', [thisFileName, '.m']);
return






densR = 2400;   VpR = 3400;   VsR = VpR*0.5;
epsR = 0;   gamR = 0;   delR = 0;   axAzimR = 0;   axTiltR = 0;
densT = 2500;   VpTv = 4152.6;   VsTv = 1793.8;
epsTv = -0.1;   delTv = -0.1;   gamTv = 0.1;   axAzimT = 0;   axTiltT = 90;
epsT = -epsTv/(1+2*epsTv);  gamT = gamTv;
VpT = VpTv/sqrt(1+2*epsT); VsT = VsTv/sqrt(1+2*gamT);
f = 1-VsT^2/VpT^2;
delT = delTv*(1+2*epsT)*(1+2*epsT/f) + 2*epsT*(1+epsT/f);

CijR0 = thomsen2cij([VpR, VsR, epsR, delR, gamR]);
CijT0 = thomsen2cij([VpT, VsT, epsT, delT, gamT]);

rotAngleR = u2u([axAzimR, axTiltR, 0]', 'deg', 'rad');
rotAngleT = u2u([axAzimT, axTiltT, 0]', 'deg', 'rad');
% Apply rotation  
if ~isempty(find(rotAngleR ~= 0, 1))
    % Rotate the stiffness matrix
    rotMatrix = loc2glb(rotAngleR);
    CijR = bond(CijR0, rotMatrix);
else
    CijR = CijR0;
end
if ~isempty(find(rotAngleT ~= 0, 1)) 
    % Rotate the stiffness matrix
    rotMatrix = loc2glb(rotAngleT);
    CijT = bond(CijT0, rotMatrix);
else
    CijT = CijT0;
end

% close all;
fig = figure;
set(fig,'units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
subplot('Position',[0.05 0.09 0.90 0.88]); hold on; grid on; grid minor;
colorBase = {'r','b','g','k','m','c','y',[0.5 0.5 0.5],'r','b'};
i = 1;
plotMatr = [];
for angInc = 0:10:90
    fprintf('incident angle %.0f deg\n',angInc);
    temp2angAzim = [0:3:360];
    RToutEstim = RTcoefModifForInverse(CijR,CijT,densR,densT,deg2rad(temp2angAzim),repmat(deg2rad(angInc),1,numel(temp2angAzim)),1,1);
    plot(temp2angAzim,abs(RToutEstim),'-','Color',colorBase{i},'LineWidth',2);
    i = i+1;
end
legend(strsplit(num2str(0:10:90)));
xlim([0,360]); xlabel("azimuth angles, °");
ylim([0,1]); ylabel("R_{PP}");


