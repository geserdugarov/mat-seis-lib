%% *setCij*
% Specify the stiffness matrices in model layers 

%% 
% *NB:*
%
% * This script needs to be tailored to each specific data set
% * Velocities and stiffness coefficients should obey the selected length and time units

%%
% *Input:*

% anisType  - [noLayer, 3][char] character array specifying the layer symmetry:
%             . ISO - isotropy, 
%             . VTI - transverse isotropy (the symmetry axis is not necessarily vertical), 
%             . ORT - orthotropy,
%             . MNC - monoclinic symmetry, and
%             . TRI - triclinic symmetry
% (*) All parameters are specified in the principal axes first and then rotated if necessary

%%
% *Output:* 

% CijRot    - [6, 6, noLayer] array of the stiffness matrices in global coordinates
% CijSym    - [6, 6, noLayer] array of the stiffness matrices in principal axes
% rotAngle  - [3, noLayer] array of three angles specifying the rotation
%             rotAngle = [azim, tilt, azx1], where  
%             azim - azimuth of the rotated z1-axis with respect to the crystallographic x-axis
%             tilt - tilt of the rotated z1-axis with respect to the crystallographic z-axis
%             azx1 - azimuth of the rotated x1-axis with respect to the crystallographic x-axis
% rotMatrix - [3, 3, noLayer] array of rotation matrices given by
%             rotMatrix(:,:,ilayer) = loc2glb(rotAngle(:,ilayer))  
% density   - [1, noLayer] array of the densities
% noLayer   = size(anisType, 1) - the number of model layers  

%%
% *Author:* Vladimir Grechka 2012 - 2015

%%
function [CijRot, CijSym, rotAngle, rotMatrix, density, noLayer] = setCij(anisType)
%% Settings  
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

noLayer   = size(anisType, 1);  
density   = zeros(1, noLayer);
CijRot    = zeros(6, 6, noLayer);
CijSym    = zeros(6, 6, noLayer);
rotAngle  = zeros(3, noLayer);     rotMatrix = zeros(3, 3, noLayer);

%% Input anisotropic parameters
for ilayer=1:noLayer
    % Initialize the rotation matrix
    rotMatrix(:,:,ilayer) = eye(3);

    global modelflag;
    if modelflag==1
        switch ilayer
            case 1
%                 Vp0 = 4000;   Vs0 = 2000;
%                 epsilon = 0.15;   delta =  0.05;   gamma = 0.15;
%                 rotAngle(:,ilayer) = u2u([45, 90, 0]', 'deg', 'rad');
                dens = 2500;
                Vp0 = 3000;   Vs0 = 1500;
                epsilon = 0.1;   delta =  0.05;   gamma = 0.1;
                rotAngle(:,ilayer) = u2u([0, 0, 0]', 'deg', 'rad');
            case 2
%                 Vp0 = 3000;   Vs0 = 1500;
%                 epsilon = 0.1;   delta =  0.05;   gamma = 0.1;
%                 rotAngle(:,ilayer) = u2u([45, 90, 0]', 'deg', 'rad');
                dens = 2500;
                Vp0 = 4000;   Vs0 = 2000;
                epsilon = 0.15;   delta =  0.05;   gamma = 0.15;
                rotAngle(:,ilayer) = u2u([45, 90, 0]', 'deg', 'rad');
            case 3
%                 Vp0 = 6000;   Vs0 = 3000;
%                 epsilon = 0.1;   delta =  0.05;   gamma = 0.1;
%                 rotAngle(:,ilayer) = u2u([45, 90, 0]', 'deg', 'rad');
                dens = 2500;
                Vp0 = 5000;   Vs0 = 3000;
                epsilon = 0.1;   delta =  0.05;   gamma = 0.1;
                rotAngle(:,ilayer) = u2u([-45, 90, 0]', 'deg', 'rad');
        end
    elseif modelflag==2
        switch ilayer
          case 1
            dens = 2500;
            Vp0 = 3000;   Vs0 = 1500;
            epsilon = 0.05;   delta =  0.02;   gamma = 0.05;
            rotAngle(:,ilayer) = u2u([0, 0, 0]', 'deg', 'rad');
          case 2
            dens = 2500;
            Vp0 = 6000;   Vs0 = 3000;
            epsilon = 0.1;   delta =  0.05;   gamma = 0.1;
            rotAngle(:,ilayer) = u2u([45, 90, 0]', 'deg', 'rad');
          case 3
            dens = 2500;
            Vp0 = 5000;   Vs0 = 2500;
            epsilon = 0.05;   delta =  0.02;   gamma = 0.05;
            rotAngle(:,ilayer) = u2u([0, 0, 0]', 'deg', 'rad');
        end
    elseif modelflag==3
        switch ilayer
              case 1
                dens = 2500;
                Vp0 = 3000;   Vs0 = 1500;
                epsilon = 0.0;   delta =  0.0;   gamma = 0.0;
                rotAngle(:,ilayer) = u2u([0, 0, 0]', 'deg', 'rad');
              case 2
                dens = 2500;
                Vp0 = 4000;   Vs0 = 2000;
                epsilon = 0.0;   delta =  0.0;   gamma = 0.0;
                rotAngle(:,ilayer) = u2u([0, 0, 0]', 'deg', 'rad');
              case 3
                dens = 2500;
                Vp0 = 5000;   Vs0 = 3000;
                epsilon = 0.0;   delta =  0.0;   gamma = 0.0;
                rotAngle(:,ilayer) = u2u([0, 0, 0]', 'deg', 'rad');
        end
	% model from Pavel Lykhin
    elseif modelflag==4
        switch ilayer
              case 1
                dens = 2200;
                Vp0 = 4660;   Vs0 = 2330;
                epsilon = 0.0;   delta =  0.0;   gamma = 0.0;
                rotAngle(:,ilayer) = u2u([0, 0, 0]', 'deg', 'rad');
              case 2
                dens = 2600;
                VpV = 5920; VsV = 2900;
                epsV = -0.039;  deltaV = -0.194; gamma = 0.253;
                epsilon = -epsV/(1+2*epsV);
                Vp0 = VpV/sqrt(1+2*epsilon);
                Vs0 = VsV/sqrt(1+2*gamma);
                f = 1-(Vs0/Vp0)^2;
                delta = deltaV*(1+2*epsilon)*(1+2*epsilon/f)+2*epsilon*(1+epsilon/f);
%                 Vp0 = 5684;   Vs0 = 2364;
%                 epsilon = 0.0423;   delta =  -0.143;   gamma = 0.253;
                rotAngle(:,ilayer) = u2u([60, 90, 0]', 'deg', 'rad');
              case 3
                dens = 2500;
                Vp0 = 4200;   Vs0 = 1760;
                epsilon = 0.0;   delta =  0.0;   gamma = 0.0;
                rotAngle(:,ilayer) = u2u([0, 0, 0]', 'deg', 'rad');
        end
	% Srednee Priob'e
    elseif modelflag==5
        switch ilayer
              case 1
                dens = 2000;
                Vp0 = 1590;   Vs0 = 610;
                epsilon = 0.0;   delta =  0.0;   gamma = 0.0;
                rotAngle(:,ilayer) = u2u([0, 0, 0]', 'deg', 'rad');
              case 2
                dens = 2100;
                Vp0 = 1820;   Vs0 = 610;
                epsilon = 0.0;   delta =  0.0;   gamma = 0.0;
                rotAngle(:,ilayer) = u2u([0, 0, 0]', 'deg', 'rad');
              case 3
                dens = 2200;
                Vp0 = 2120;   Vs0 = 960;
                epsilon = 0.0;   delta =  0.0;   gamma = 0.0;
                rotAngle(:,ilayer) = u2u([0, 0, 0]', 'deg', 'rad');
              case 4
                dens = 2300;
                Vp0 = 2760;   Vs0 = 1070;
                epsilon = 0.0;   delta =  0.0;   gamma = 0.0;
                rotAngle(:,ilayer) = u2u([0, 0, 0]', 'deg', 'rad');
              case 5
                dens = 2300;
                Vp0 = 3260;   Vs0 = 1360;
                epsilon = 0.0;   delta =  0.0;   gamma = 0.0;
                rotAngle(:,ilayer) = u2u([0, 0, 0]', 'deg', 'rad');
              case 6
                dens = 2400;
                Vp0 = 3740;   Vs0 = 1700;
                epsilon = 0.0;   delta =  0.0;   gamma = 0.0;
                rotAngle(:,ilayer) = u2u([0, 0, 0]', 'deg', 'rad');
              case 7
                dens = 2300;
                Vp0 = 3280;   Vs0 = 1810;
                epsilon = 0.0;   delta =  0.0;   gamma = 0.0;
                rotAngle(:,ilayer) = u2u([0, 0, 0]', 'deg', 'rad');
              case 8
                dens = 2400;
                VpV = 3600; VsV = 1900;
                epsV = -0.1;  deltaV = -0.06; gamma = 0.1;
                epsilon = -epsV/(1+2*epsV);
                Vp0 = VpV/sqrt(1+2*epsilon);
                Vs0 = VsV/sqrt(1+2*gamma);
                f = 1-(Vs0/Vp0)^2;
                delta = deltaV*(1+2*epsilon)*(1+2*epsilon/f)+2*epsilon*(1+epsilon/f);
                rotAngle(:,ilayer) = u2u([75, 90, 0]', 'deg', 'rad');
        end
	% Verhnechonskoe
    elseif modelflag==6
        switch ilayer
              case 1
                dens = 2100;
                Vp0 = 2260;   Vs0 = 1020;
                epsilon = 0.0;   delta =  0.0;   gamma = 0.0;
                rotAngle(:,ilayer) = u2u([0, 0, 0]', 'deg', 'rad');
              case 2
                dens = 2700;
                Vp0 = 6500;   Vs0 = 3380;
                epsilon = 0.0;   delta =  0.0;   gamma = 0.0;
                rotAngle(:,ilayer) = u2u([0, 0, 0]', 'deg', 'rad');
              case 3
                dens = 2600;
                Vp0 = 5680;   Vs0 = 2850;
                epsilon = 0.0;   delta =  0.0;   gamma = 0.0;
                rotAngle(:,ilayer) = u2u([0, 0, 0]', 'deg', 'rad');
              case 4
                dens = 2600;
                Vp0 = 6000;   Vs0 = 2700;
                epsilon = 0.0;   delta =  0.0;   gamma = 0.0;
                rotAngle(:,ilayer) = u2u([0, 0, 0]', 'deg', 'rad');
              case 5
                dens = 2600;
                Vp0 = 4880;   Vs0 = 2540;
                epsilon = 0.0;   delta =  0.0;   gamma = 0.0;
                rotAngle(:,ilayer) = u2u([0, 0, 0]', 'deg', 'rad');
              case 6
                dens = 2600;
                Vp0 = 5740;   Vs0 = 2300;
                epsilon = 0.0;   delta =  0.0;   gamma = 0.0;
                rotAngle(:,ilayer) = u2u([0, 0, 0]', 'deg', 'rad');
              case 7
                dens = 2200;
                Vp0 = 4660;   Vs0 = 2330;
                epsilon = 0.0;   delta =  0.0;   gamma = 0.0;
                rotAngle(:,ilayer) = u2u([0, 0, 0]', 'deg', 'rad');
              case 8
                dens = 2600;
                VpV = 5920; VsV = 2900;
                epsV = -0.1;  deltaV = -0.06; gamma = 0.1;
                epsilon = -epsV/(1+2*epsV);
                Vp0 = VpV/sqrt(1+2*epsilon);
                Vs0 = VsV/sqrt(1+2*gamma);
                f = 1-(Vs0/Vp0)^2;
                delta = deltaV*(1+2*epsilon)*(1+2*epsilon/f)+2*epsilon*(1+epsilon/f);
                rotAngle(:,ilayer) = u2u([60, 90, 0]', 'deg', 'rad');
              case 9
                dens = 2500;
                Vp0 = 4180;   Vs0 = 1760;
                epsilon = 0.0;   delta =  0.0;   gamma = 0.0;
                rotAngle(:,ilayer) = u2u([0, 0, 0]', 'deg', 'rad');
              case 10
                dens = 2700;
                Vp0 = 6200;   Vs0 = 3530;
                epsilon = 0.0;   delta =  0.0;   gamma = 0.0;
                rotAngle(:,ilayer) = u2u([0, 0, 0]', 'deg', 'rad');
        end
	% Turia 1
    elseif modelflag==7
        switch ilayer
              case 1
                dens = 1980;
                Vp0 = 1750;   Vs0 = 880;
                epsilon = 0.0;   delta =  0.0;   gamma = 0.0;
                rotAngle(:,ilayer) = u2u([0, 0, 0]', 'deg', 'rad');
              case 2
                dens = 2050;
                Vp0 = 2100;   Vs0 = 1000;
                epsilon = 0.0;   delta =  0.0;   gamma = 0.0;
                rotAngle(:,ilayer) = u2u([0, 0, 0]', 'deg', 'rad');
              case 3
                dens = 2330;
                Vp0 = 2900;   Vs0 = 1500;
                epsilon = 0.0;   delta =  0.0;   gamma = 0.0;
                rotAngle(:,ilayer) = u2u([0, 0, 0]', 'deg', 'rad');
              case 4
                dens = 2550;
                VpV = 4200; VsV = 2400;
                epsV = -0.01;  deltaV = 0.08; gamma = 0.01;
                epsilon = -epsV/(1+2*epsV);
                Vp0 = VpV/sqrt(1+2*epsilon);
                Vs0 = VsV/sqrt(1+2*gamma);
                f = 1-(Vs0/Vp0)^2;
                delta = deltaV*(1+2*epsilon)*(1+2*epsilon/f)+2*epsilon*(1+epsilon/f);
                rotAngle(:,ilayer) = u2u([30-90, 90, 0]', 'deg', 'rad');
%               case 4
%                 dens = 2550;
%                 Vp0 = 4200; Vs0 = 2400;
%                 epsilon = 0.0;   delta =  0.0;   gamma = 0.0;
%                 rotAngle(:,ilayer) = u2u([0, 0, 0]', 'deg', 'rad');
%               case 5
%                 dens = 2730;
%                 VpV = 4770; VsV = 3000;
%                 epsV = -0.01;  deltaV = 0.05; gamma = 0.01;
%                 epsilon = -epsV/(1+2*epsV);
%                 Vp0 = VpV/sqrt(1+2*epsilon);
%                 Vs0 = VsV/sqrt(1+2*gamma);
%                 f = 1-(Vs0/Vp0)^2;
%                 delta = deltaV*(1+2*epsilon)*(1+2*epsilon/f)+2*epsilon*(1+epsilon/f);
%                 rotAngle(:,ilayer) = u2u([150-90, 90, 0]', 'deg', 'rad');
              case 5
                dens = 2730;
                Vp0 = 4770; Vs0 = 3000;
                epsilon = 0.0;   delta =  0.0;   gamma = 0.0;
                rotAngle(:,ilayer) = u2u([0, 0, 0]', 'deg', 'rad');
        end
    else
        switch ilayer
              case 1
                dens = 2500;
                Vp0 = 2500;   Vs0 = 1300;  
                epsilon = 0.150;   delta =  0.050;   gamma = -0.1;
                rotAngle(:,ilayer) = u2u([10, 10, 0]', 'deg', 'rad'); 
              case 2
                 dens = 2500;
                 Vp0 = 3000;   Vs0 = 1700;
                 epsilon = 0.200;   delta =  0.10;   gamma = -0.1;
                 rotAngle(:,ilayer) = u2u([20, 10, 0]', 'deg', 'rad');
              case 3
                dens = 2500;
                Vp0 = 3500;   Vs0 = 2000;
                epsilon = 0.250;   delta =  0.100;   gamma = -0.2;
                rotAngle(:,ilayer) = u2u([30, 10, 0]', 'deg', 'rad');
              case 4
                dens = 2500;
                Vp0 = 3800;   Vs0 = 2200;
                epsilon = 0.050;   delta =  0.01;   gamma = -0.2;
                rotAngle(:,ilayer) = u2u([40, 10, 0]', 'deg', 'rad');
        end
    end
    density(ilayer) = dens;
    
    epsilon1 = epsilon;   epsilon2 = epsilon+0.05;                  % ORT setups
    delta1 = delta;       delta2 = delta+0.05;    delta3 = -0.02; 
    gamma1 = gamma;       gamma2 = gamma+0.05;

    %%
    if strcmp(anisType(ilayer,:), 'ISO') == 1
        % Input Vp and Vs
        C0 = thomsen2cij([Vp0, Vs0, 0, 0, 0]);
    end
    
    %%
    if strcmp(anisType(ilayer,:), 'VTI') == 1
        % Input Thomsen's parameters 
        C0 = thomsen2cij([Vp0, Vs0, epsilon, delta, gamma]);
    end

    %%
    if strcmp(anisType(ilayer,:), 'ORT') == 1
        % Input Tsvankin's parameters 
        C0 = tsvankin2cij([Vp0, Vs0, epsilon1, epsilon2, delta1, delta2, delta3, ...
                    gamma1, gamma2]);  
        rotAngle(:,ilayer) = u2u([30, 30, 50]', 'deg', 'rad');
    end

    %%
    if strcmp(anisType(ilayer,:), 'MNC') == 1
        % Input Grechka's parameters
        zeta1 = 0.01;   zeta2 = 0.02;   zeta3 = 0.03;
        C0 = grechka2cij([Vp0, Vs0, epsilon1, epsilon2, delta1, delta2, delta3, ...
                      gamma1, gamma2, zeta1, zeta2, zeta3]);
        rotAngle(:,ilayer) = u2u([60, 30, 80]', 'deg', 'rad');
    end

    %%
    if strcmp(anisType(ilayer,:), 'TRI') == 1
        C1 = 0.1*[107, 18, 19, 12, 8, -1, 111, 38, -2, -3, 7, 101, 8, -7, -6, 27, 1, 4, 42, -16, 51]; 
        C0 = vec2mat(C1);
    end

    %% Check the stability conditions
    [flag, ~] = isStable(C0, 0); 
    if flag == 0
        fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
        fprintf('>>> The stability conditions are violated in layer %g \n \n', ilayer);
          error('>>> STOP');
    end  

    %% Write out the stiffness tensor in crystallographic coordinate frame
    CijSym(:,:,ilayer) = C0; 
    
    % Apply rotation  
    if isempty(find(rotAngle(:,ilayer) ~= 0, 1)) ~= 1 
        % Rotate the stiffness matrix
        rotMatrix(:,:,ilayer) = loc2glb(rotAngle(:,ilayer));
        CijRot(:,:,ilayer) = bond(C0, rotMatrix(:,:,ilayer));
    else
        CijRot(:,:,ilayer) = C0;
    end

end  % of loop over ilayer

end    % of the function