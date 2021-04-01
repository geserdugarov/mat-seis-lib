%% *art3D*
% Two-point ray tracing in 3D anisotropic media 

%%
% *Input:*

% model.         - structure containing the following arrays:
%     interface  - [5, noInt] interface array specified by function 'setInterface'
%          noInt = size(interface,2) - the number of interfaces
%     fault      - [5, noFlt] fault array specified by function 'setFault'
%          noFlt = size(fault,2) - the number of faults
%     Cij.global - [6, 6, noInt+1] array of the stiffness matrices of layers (in the global
%                  coordinate frame) ordered in increasing depth
%     waveCode   - [4, noWave] wave-code array specified by function 'setWaveCode'
%         noWave - the number of wave codes
%     xSou       - [3, noSou] array of the source coordinates
%          noSou - the number of sources
%     xRec       - [3, noRec] array of the receiver coordinates
%          noRec - the number of receivers
%
% artFlags.      - structure containing the following quantities:
%     flagWAA    - flag indicating whether calculations should be performed for weak
%                  (flagWAA = 'WA') or strong (flagWAA = 'SA') anisotropy  
%     flagDT     - [scalar] flag indicating whether to calculate the Frechet derivatives of 
%                  traveltimes (flagDT = 1) or skip the calculation (flagDT = 0)
%     flagDU     - [scalar] flag indicating whether to calculate the Frechet derivatives of 
%                  polarization vectors (flagDU = 1) or skip the calculation (flagDU = 0)

%%
% *Output:*

% art.           - structure containing the following fields:
%     time       - [noTime, 1] array of computed traveltimes 
%         noTime = noSou*noRec*noWave 
%     rcode      - [4, 2*noInt+2+noFlt, noTime] array of ray codes created by function 
%                  'produceRayCode'
%     tseg       - [2*noInt+1+noFlt, noTime] array of propagation times along the ray segments
%     traj       - [3, 2*noInt+2+noFlt, noTime] array of coordinates of points along a ray 
%                  trajectory sorted as [source, intersections with interfaces and faults, receiver] 
%                  (*) The second dimension '2*noInt+2+noFlt' implies that reflections from 
%                      interfaces and faults are allowed but multiples are not  
%     n          - [3, 2*noInt+1+noFlt, noTime] array of the unit wavefront normal vectors along   
%                  the ray segments
%     r          - [3, 2*noInt+1+noFlt, noTime] array of the unit ray vectors along the ray segments
%     p          - [3, 2*noInt+1+noFlt, noTime] array of the slowness vectors along the ray segments
%                  (*) Note: p ~= n/Vph in the WAA
%     U          - [3, 2*noInt+1+noFlt, noTime] array of the unit polarization vectors along 
%                  the ray segments
%     Vph        - [2*noInt+1+noFlt, noTime] array of phase velocities along the ray segments
%     Vgr        - [2*noInt+1+noFlt, noTime] array of group velocities along the ray segments
%     Kp         - [2*noInt+1+noFlt, noTime] array of the Gaussian curvatures of the slowness 
%                  surfaces along the ray segments
%     dtdm       - [noTime, 21*(noInt+1) + 3*noSou + 3*noRec + 3*noInt + noSou + 3] array of 
%                  the Frechet derivatives of artTime with respect to
%                  . 21 stiffness components,
%                  . 3 source coordinates,
%                  . 3 receiver coordinates,  
%                  . 3 interface parameters [depth and the x- and y-components of its normal], 
%                    the event-origin times, and
%                  . 3 fault parameters             
%     dUdm       - [3*noTime, 21*(noInt+1) + 3*noSou + 3*noRec + 3*noInt + 3*noFlt + noSou] array 
%                  of the Frechet derivatives of artU(:,end,:) with respect to the model parameters 
%                  listed above

%%
% *Author:* Vladimir Grechka 2012 - 2014

%% 
% *Known issues:* 
%
% * fault reflections might be improperly handled; more testing is required  
% * to be implemented: 
%   - derivatives |art.dUdm| 
%   - check of Snell's law at the model interfaces ... how much violation is acceptable?

%%
function [art] = art3D(model, artFlags, timePicks)
%% Settings 
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

%% Unpack the input structures into local variables
interface = model.interface;    fault  = model.fault;      Cij    = model.Cij.global;
waveCode  = model.waveCode;     xSou   = model.xSou;       xRec   = model.xRec;
flagWAA   = artFlags.flagWAA;   flagDT = artFlags.flagDT;  flagDU = artFlags.flagDU;

%% Check the calculation mode 
if strcmp(flagWAA, 'SA') == 0 && strcmp(flagWAA, 'WA') == 0
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Flag flagWAA = ''%s'' is unsupported \n \n', flagWAA); 
    fprintf('    Supported flags are flagWAA = ''SA'' and flagWAA = ''WA'' \n'); 
      error('>>> STOP');
end

%% Check the consistence of arrays |interface| and |Cij|
if size(Cij, 3) ~= size(interface, 2) + 1  
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf(['>>> The number of layers (%g) is not equal to the number of interfaces ', ...
             '(%g) + 1 \n'], [size(Cij, 3), size(interface, 2)]); 
    fprintf('    Please check the functions ''setCij'' and ''setInterface'' \n \n'); 
      error('>>> STOP');
end

%% Optimization options
tolFun = 1.e-8;  tolX = 1.e-6;  tolSnell = 1.e-2;
%tolFaultOptim = 1.e-3;  tolGradCheck = 1.e-2;  
noIter = 1000;
OPTIONSunc = ...  
    optimset('LargeScale', 'on', 'Display','off', 'GradObj', 'on', ...  % 'Hessian', 'on', ...
             'DerivativeCheck', 'off', ...
             'MaxIter', noIter, 'MaxFunEvals', noIter, 'TolX', tolX, 'TolFun', tolFun);     
OPTIONScon = ...
    optimset('Algorithm', 'active-set', 'Display', 'off', 'GradObj', 'on', 'GradConstr', 'off', ...
             'DerivativeCheck', 'off', ...
             'MaxIter', noIter, 'MaxFunEvals', noIter, 'TolX', tolX, 'TolFun', tolFun);     

%% Sizes of the output arrays
noSou   = size(xSou, 2);        noRec   = size(xRec, 2);    noWave  = size(waveCode, 2);
noInt   = size(interface, 2);   noFlt   = size(fault, 2);   noTime  = noSou*noRec*noWave;
noParam = 21*(noInt+1) + 3*noSou + 3*noRec + 3*noInt + 3*noFlt + noSou;  

%% Initialize arrays
noIntFlt1 = 2*noInt + 1 + noFlt;            noIntFlt2 = noIntFlt1 + 1; 
art.time  = zeros(noTime, 1);               art.traj  = NaN(3, noIntFlt2, noTime);
art.rcode = zeros(4, noIntFlt2, noTime);    art.tseg  = zeros(noIntFlt1, noTime);       
art.n     = zeros(3, noIntFlt1, noTime);    art.Vph   = zeros(noIntFlt1, noTime);
art.U     = zeros(3, noIntFlt1, noTime);    art.Vgr   = zeros(noIntFlt1, noTime);
art.r     = zeros(3, noIntFlt1, noTime);    art.Kp    = zeros(noIntFlt1, noTime);
art.p     = zeros(3, noIntFlt1, noTime);
if flagDT == 1
    art.dtdm = zeros(noTime, noParam);  
else
    art.dtdm = [];
end
art.dUdm  = [];  % replacement for 'art.dUdm  = zeros(3*noTime, noParam);' to save memory 

artTseg1  = zeros(   noIntFlt1, noWave, noRec, noSou);
artN1     = zeros(3, noIntFlt1, noWave, noRec, noSou);
artU1     = zeros(3, noIntFlt1, noWave, noRec, noSou);
artR1     = zeros(3, noIntFlt1, noWave, noRec, noSou);
artP1     = zeros(3, noIntFlt1, noWave, noRec, noSou);
artVph1   = zeros(   noIntFlt1, noWave, noRec, noSou);
artVgr1   = zeros(   noIntFlt1, noWave, noRec, noSou);
artKp1    = zeros(   noIntFlt1, noWave, noRec, noSou);

%artTime1  = zeros(noWave, noRec, noSou);
%artTraj1  =   NaN(3, noIntFlt2, noWave, noRec, noSou);
%artRcode1 =   NaN(4, noIntFlt2, noWave, noRec, noSou);
%artDT1    = zeros(noWave, noRec, noSou, noParam);

%% The z-coefficients and normals of the interfaces and faults
[zInt, nrmInt] = planeEquation(interface);
[zFlt, nrmFlt] = planeEquation(fault);

%% Ray tracing
%parfor isou = 1:noSou
% parfor irec = 1:noRec
for irec = 1:noRec
    for isou = 1:noSou
        for iwave = 1:noWave   
%            [isou, irec, iwave]
            %% Initialize temporary arrays for parallel computations
            flagOut = 1;
            artTraj0 = NaN(3, noIntFlt2);    artRcode0 = NaN(4, noIntFlt2);
            artN0    = zeros(3, noIntFlt1);  artR0   = artN0;    artU0  = artN0;   artP0  = artN0; 
            artVph0  = zeros(noIntFlt1, 1);  artVgr0 = artVph0;  artKp0 = artVph0;
            artTseg0 = artVph0;
            artDT0   = zeros(1, noParam);           
%            artDU0   = zeros(3, noParam);
            
            % Check whether timePicks is valid
            itime = noRec*noWave*(isou - 1) + noWave*(irec - 1) + iwave; 
            if isempty(timePicks) == 1  ||  isnan(timePicks(itime)) == 0      
                %% Loop over potential fault reflections
                for ilayer = 1:(noInt+1)            
                    fltLayer = ilayer;
                    % Produce a ray code
                    [rayCode, sou3, rec3, zIntCur, nrmIntCur, CijCur] = ...
                        produceRayCode(xSou(:,isou), xRec(:,irec), zInt, zFlt, nrmInt, nrmFlt, ...
                                       Cij, waveCode(:,iwave), fltLayer);
                    souCur = [xSou(1:2,isou); sou3];    % the current source and receiver coordinates
                    recCur = [xRec(1:2,irec); rec3];
                    noSeg = size(rayCode, 2);           % the number of segments of the ray trajectory
                    rayType = rayCode(1,:)';            % wave types along the trajectory

                    % Initialize the ray trajectory
                    [xy0] = initRayTraj(souCur, recCur, waveCode(3,iwave), zInt, ...
                                                        waveCode(4,iwave), zFlt, fltLayer, rayCode);
                    if noSeg == 1
                        % Single ray segment  
                        [artTimeCur, ~] = ...
                            rayTracingOptim(xy0, rayType, souCur, recCur, zIntCur, CijCur, flagWAA);
                        xy = xy0;
                        break;  % out of the fault segment loop
                            
                    else
                        % Optimization for the ray trajectory
                        if sum(rayCode(4,:)) == 0                       % unconstrained optimization                             
                            [xy, artTimeCur] = ...
                                fminunc('rayTracingOptim', xy0, ...
                                    OPTIONSunc, rayType, souCur, recCur, zIntCur, CijCur, flagWAA);
                             break;  % exit the fault-segment loop

                        else                                            % constrained optimization
                            [Acon, Bcon] = faultConstraint(xy0, zInt, zFlt, fltLayer, rayCode);
                            [xy, artTimeCur] = ...
                                fmincon('rayTracingOptim', xy0, Acon, Bcon, [], [], [], [], [], ...
                                    OPTIONScon, rayType, souCur, recCur, zIntCur, CijCur, flagWAA); 
                             break;  % exit the fault-segment loop
                        end
                    end
                end    % end of loop over the fault segments

                %% Fill out output arrays    
                if flagOut ~= 0
                    % Get trajectory information
                    [xyz, tOut, nOut, rOut, pOut, UOut, VOut, dVOut, gOut, dgOut, KpOut, ...
                        rayTypeOut, flagOut] = ...
                        rayAttributes(xy, rayType, souCur, recCur, zIntCur, CijCur, flagWAA);
                     rayCode(1,:) = rayTypeOut'; 
                     [~, flagSnell] = isSnell(pOut, nrmIntCur, tolSnell);
                     if flagSnell == 0
                         fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, ...
                                 [thisFileName, '.m']));
                         fprintf(['>>> Violation of Snell''s law ', ...
                                  'for isou = %g,  irec = %g,  iwave = %g \n'], isou, irec, iwave);
                        % Set the error flag
                        artTimeCur = NaN;  flagOut = 0;                     
                     end
                end
            else
                % timePicks(itime) = NaN
                flagOut = 0;
            end
                  
            if flagOut == 0
                artTimeCur = NaN;          
                if flagDT == 1        
                    artDT1(iwave, irec, isou, :) = zeros(1, noParam);
                end
                artTime1(iwave, irec, isou)        = artTimeCur;
                artRcode1(:, :, iwave, irec, isou) = artRcode0;
                artN1(:, :, iwave, irec, isou)     = artN0;
                artU1(:, :, iwave, irec, isou)     = artU0;
                artR1(:, :, iwave, irec, isou)     = artR0;
                artP1(:, :, iwave, irec, isou)     = artP0;
                artTseg1(:, iwave, irec, isou)     = artTseg0;
                artTraj1(:, :, iwave, irec, isou)  = artTraj0;
                artVph1(:, iwave, irec, isou)      = artVph0;
                artVgr1(:, iwave, irec, isou)      = artVgr0;
                artKp1(:, iwave, irec, isou)       = artKp0;
            else
                % Fill in output arrays
                artTime1(iwave, irec, isou)        = artTimeCur;
                artRcode1(:, :, iwave, irec, isou) = copyMat(rayCode, artRcode0);
                artN1(:, :, iwave, irec, isou)     = copyMat(nOut,    artN0);
                artU1(:, :, iwave, irec, isou)     = copyMat(UOut,    artU0);
                artR1(:, :, iwave, irec, isou)     = copyMat(rOut,    artR0);
                artP1(:, :, iwave, irec, isou)     = copyMat(pOut,    artP0);
                artTseg1(:, iwave, irec, isou)     = copyMat(tOut',   artTseg0);
                artTraj1(:, :, iwave, irec, isou)  = copyMat(xyz',    artTraj0);
                artVph1(:, iwave, irec, isou)      = copyMat(VOut',   artVph0);
                artVgr1(:, iwave, irec, isou)      = copyMat(gOut',   artVgr0);
                artKp1(:, iwave, irec, isou)       = copyMat(KpOut',  artKp0);          
            end
                
            %% Frechet derivatives of the traveltime
            if flagDT == 1 && flagOut ~= 0
                % dtime/dCij
                for iseg = 1:noSeg
                    if strcmp(flagWAA, 'SA') == 1
                        [~, dgVec] = dGdCij(CijCur(:,:,iseg), nOut(:,iseg)/VOut(iseg), ...
                                   UOut(:,iseg), gOut(:,iseg));
                    else
                        [~, dgVec] = dVdCij(nOut(:,iseg), VOut(iseg), UOut(:,iseg));
                    end
                    ilayer = 21*rayCode(2,iseg);
                    artDT0(1, ilayer-20:ilayer) = artDT0(1, ilayer-20:ilayer) - ...
                            norm(xyz(iseg+1,:) - xyz(iseg,:))*dgVec/gOut(iseg)^2;
                end

                % dtime/dsou
                count = 21*(noInt+1);
                artDT0(1, count + 3*isou - 2 : count + 3*isou) = -pOut(:,1);

                % dtime/drec
                count = 21*(noInt+1) + 3*noSou;
                artDT0(1, count + 3*irec - 2 : count + 3*irec) = pOut(:,noSeg);
                
                % dtime/dint and dtime/dflt
                count = 21*(noInt+1) + 3*noSou + 3*noRec;

                if noSeg > 1
                    for iseg = 1:(noSeg-1)
                        if rayCode(4,iseg) ~= 0
                            n1 = nrmFlt(:,rayCode(4,iseg));
                            step1 = noInt + rayCode(4,iseg);
                            dn = [0, 0, 0; 1, 0, -n1(1)/n1(3); 0, 1, -n1(2)/n1(3)];
                            factor1 = ([n1(3,1); 0; 0] + dn*(fault(1:3, step1) - xyz(iseg+1,:)'))';
                        else
                            n1 = nrmInt(:,rayCode(3,iseg));
                            step1 = rayCode(3,iseg);
                            dn = [0, 0, 0; 1, 0, -n1(1)/n1(3); 0, 1, -n1(2)/n1(3)];
                            factor1 = ([n1(3,1); 0; 0] + dn*(interface(1:3, step1) - xyz(iseg+1,:)'))';
                        end                    
                        artDT0(1, count + 3*step1 - 2 : count + 3*step1) = ...
                        artDT0(1, count + 3*step1 - 2 : count + 3*step1) + ...
                            n1'*(dTdR((xyz(iseg+1,:) - xyz(iseg,:))', ...
                                    VOut(iseg),   dVOut(:,iseg),   dgOut(:,iseg),   flagWAA) - ...
                                 dTdR((xyz(iseg+2,:) - xyz(iseg+1,:))', ...
                                    VOut(iseg+1), dVOut(:,iseg+1), dgOut(:,iseg+1), flagWAA))* ...
                            factor1;
                    end
                elseif noInt > 0 || noFlt > 0  % ray intersects neither interfaces nor faults ->
                    artDT0(1, count + 1 : count + 3*noInt + 3*noFlt) = 0;        % fill in zeros 
                else
                    % Do nothing. The model has no interfaces and no faults.
                end
                    
                % dtime/dtau
                count = 21*(noInt+1) + 3*noSou + 3*noRec + 3*noInt + 3*noFlt;
                artDT0(1, count + isou) = 1;
                artDT1(iwave, irec, isou, :) = artDT0;
            end    % of computation of the Frechet derivatives of traveltimes

            %% Frechet derivatives of the polarization
            if flagDU == 1 && flagOut ~= 0
                % Do nothing. This part is to be completed when the Frechet derivatives of artU 
                % are available.
            end

        end    % of loop over iwave
    end    % of loop over irec
end    % of loop over isou

%% Fill in arrays of the output structure
for isou = 1:noSou
    for irec = 1:noRec
        for iwave = 1:noWave
            itime = noRec*noWave*(isou - 1) + noWave*(irec - 1) + iwave; 
            art.time(itime, 1)     = artTime1(iwave, irec, isou);
            art.rcode(:, :, itime) = artRcode1(:, :, iwave, irec, isou);
            art.tseg(:, itime)     = artTseg1(:, iwave, irec, isou);
            art.traj(:, :, itime)  = artTraj1(:, :, iwave, irec, isou);
            art.n(:, :, itime)     = artN1(:, :, iwave, irec, isou);
            art.r(:, :, itime)     = artR1(:, :, iwave, irec, isou);
            art.p(:, :, itime)     = artP1(:, :, iwave, irec, isou);
            art.U(:, :, itime)     = artU1(:, :, iwave, irec, isou);
            art.Vph(:, itime)      = artVph1(:, iwave, irec, isou);
            art.Vgr(:, itime)      = artVgr1(:, iwave, irec, isou);
            art.Kp(:, itime)       = artKp1(:, iwave, irec, isou);
            if flagDT == 1  
                art.dtdm(itime, :)  = artDT1(iwave, irec, isou, :);
            end
        end
    end
end        
            
end    % of the function

%%
