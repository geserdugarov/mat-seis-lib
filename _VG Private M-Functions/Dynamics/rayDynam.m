%% *rayDynam*
% Ray dynamic calculation after kynematics were done 

%%
% *Input:*

% model.         - structure containing the following arrays:
%     interface  - [5, noInt] interface array specified by function 'setInterface'
%          noInt = size(interface,2) - the number of interfaces
%     fault      - [5, noFlt] fault array specified by function 'setFault'
%          noFlt = size(fault,2) - the number of faults
%     Cij.global - [6, 6, noInt+1] array of the stiffness matrices of layers (in the global
%                  coordinate frame) ordered in increasing depth
%     density    - [1, noInt+1] array of the densities
%     waveCode   - [4, noWave] wave-code array specified by function 'setWaveCode'
%         noWave - the number of wave codes
%     xSou       - [3, noSou] array of the source coordinates
%          noSou - the number of sources
%     xRec       - [3, noRec] array of the receiver coordinates
%          noRec - the number of receivers

% artOut.        - structure containing the following arrays:
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
% *Output:*

% dynamOut.      - structure containing the following fields:
%     RTcoef     - [noRT, noRays] reflection/refraction coefficients
%     dpdA       - [3, noRaySeg, noRays] for geometrical spreading
%     dpdB       - [3, noRaySeg, noRays] for geometrical spreading  
%     dxdA       - [3, noRaySeg, noRays] for geometrical spreading
%     dxdB       - [3, noRaySeg, noRays] for geometrical spreading
%     L          - [1, noRays] geometrical spreading
%     P          - [1, noRays] adding reflection/refraction
%     Amp        - [1, noRays] amplitudes in the receivers

%%
% *Author:* Geser Dugarov 2017

%% 
% *Known issues:* 
%
% * the script was not been tested
% * only for planar horizontal layers

%%
function [dynam] = rayDynam(model, artOut)
%% Settings 
[~, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

[~, ~, ~, t2m] = indexesVMT(4); % for Cij from [6,6] to [3,3,3,3]

%% Unpack the input structures into local variables
rcode = artOut.rcode;
rho   = model.density;
Aij   = model.Cij.global;
Cij   = zeros(size(Aij));
for i = 1:numel(rho)
    Cij(:,:,i) = rho(i)*Aij(:,:,i);
end

%% Initialize arrays
noRays = size(rcode,3);
temp = rcode(1,:,1);
temp(isnan(temp)) = [];
% all rays have the same segment numbers
noRaySeg = numel(temp); % number of ray segments
clear temp;
noRT = noRaySeg - 1;

dynam.RTcoef = ones(noRT, noRays);
dynam.dpdA   = zeros(3, noRaySeg*2, noRays);
dynam.dpdB   = zeros(3, noRaySeg*2, noRays);
dynam.dxdA   = zeros(3, noRaySeg*2, noRays);
dynam.dxdB   = zeros(3, noRaySeg*2, noRays);
dynam.L      = zeros(1, noRays);
dynam.P      = zeros(1, noRays);
dynam.Amp    = zeros(1, noRays);

%% Dynamic ray tracing
% p3 for searching all R/T waves
% (R/T)r - reflection/transmission coefficient of the r-wave
syms p3 RP RS1 RS2 TP TS1 TS2;
% parfor raynum = 1:noRays
for raynum = 1:noRays
    fprintf('Ray ''%d'' from ''%d''.\n', raynum, noRays);
    raydata = rcode(:,1:noRaySeg,raynum);
    t   = artOut.tseg(1:noRaySeg,raynum);
    Vph = artOut.Vph(1:noRaySeg,raynum);
    Vgr = artOut.Vgr(1:noRaySeg,raynum);
    n   = artOut.n(:,1:noRaySeg,raynum);
    p   = artOut.p(:,1:noRaySeg,raynum);
    r   = artOut.r(:,1:noRaySeg,raynum);
    U   = artOut.U(:,1:noRaySeg,raynum);
    
    %% geometrical spreading
    % starting values
    n1 = n(:,1);
    dndA = [n1(1)*n1(3)/sqrt(1-n1(3)^2); n1(2)*n1(3)/sqrt(1-n1(3)^2); -sqrt(1-n1(3)^2)];
    dndB = [-n1(2); n1(1); 0];
    dynam.dpdA(:,1,raynum) = 1/Vph(1).*dndA - n1./Vph(1)*dot(r(:,1)./Vgr(1),dndA); % + or - ???
    dynam.dpdB(:,1,raynum) = 1/Vph(1).*dndB - n1./Vph(1)*dot(r(:,1)./Vgr(1),dndA); % + or - ???
    dynam.dxdA(:,1,raynum) = 0;
    dynam.dxdB(:,1,raynum) = 0;
    for elem = 2:noRaySeg*2
        gk = [0; 0; 1]; % normal to the bound
        segNum = ceil(elem/2);
        layer  = raydata(2,segNum);
        wavetype = raydata(1,segNum);
        [X,V] = velPhaseU(Aij(:,:,layer),n(:,segNum));
        [~,I] = sort(X,'descend');
        if mod(elem,2) == 0
            % through layer
            dynam.dpdA(:,elem,raynum) = dynam.dpdA(:,elem-1,raynum);
            dynam.dpdB(:,elem,raynum) = dynam.dpdB(:,elem-1,raynum);
            dUdA = zeros(3,1);
            dUdB = zeros(3,1);
            BrqA = 0;
            BrqB = 0;
            for i = 1:3
                % wavetype: 1-qP, 2-qS1, 3-qS2; i order: 1-qP, 2-qSV, 1-SH
                if i ~= wavetype
                    for j = 1:3
                        for k = 1:3
                            for l = 1:3
                                for m = 1:3
                                    BrqA = BrqA + (Aij(t2m(j,k),t2m(l,m),layer) + Aij(t2m(j,l),t2m(k,m),layer))*...
                                        dynam.dpdA(k,elem,raynum)*p(l,segNum)*U(m,segNum)*V(j,I(i));
                                    BrqB = BrqB + (Aij(t2m(j,k),t2m(l,m),layer) + Aij(t2m(j,l),t2m(k,m),layer))*...
                                        dynam.dpdB(k,elem,raynum)*p(l,segNum)*U(m,segNum)*V(j,I(i));
                                end
                            end
                        end
                    end
                    BrqA = BrqA*Vph(segNum)^2/(Vph(segNum)^2-X(I(i)));
                    BrqB = BrqB*Vph(segNum)^2/(Vph(segNum)^2-X(I(i)));
                    dUdA = dUdA + BrqA*V(:,I(i));
                    dUdB = dUdB + BrqB*V(:,I(i));
                end
            end
            dvdA = zeros(3,1);
            dvdB = zeros(3,1);
            for j = 1:3
                for k = 1:3
                    for l = 1:3
                        for m = 1:3
                            dvdA(j) = dvdA(j) + ...
                               (Aij(t2m(j,k),t2m(l,m),layer) + Aij(t2m(j,m),t2m(l,k),layer))*p(l,segNum)*dUdA(k)*U(m,segNum) + ...
                                Aij(t2m(j,k),t2m(l,m),layer)*dynam.dpdA(l,elem,raynum)*U(k,segNum)*U(m,segNum);
                            dvdB(j) = dvdB(j) + ....
                               (Aij(t2m(j,k),t2m(l,m),layer) + Aij(t2m(j,m),t2m(l,k),layer))*p(l,segNum)*dUdB(k)*U(m,segNum) + ...
                                Aij(t2m(j,k),t2m(l,m),layer)*dynam.dpdB(l,elem,raynum)*U(k,segNum)*U(m,segNum);
                        end
                    end
                end
            end
            dynam.dxdA(:,elem,raynum) = dynam.dxdA(:,elem-1,raynum) + dvdA*t(segNum);
            dynam.dxdB(:,elem,raynum) = dynam.dxdB(:,elem-1,raynum) + dvdB*t(segNum);
        else
            % through bound
            dynam.dpdA(:, elem, raynum) = dynam.dpdA(:,elem-1,raynum) - ...
                gk*dot(dynam.dpdA(:,elem-1,raynum), r(:,segNum)*Vgr(segNum))/dot(gk,r(:,segNum)*Vgr(segNum));
            dynam.dpdB(:, elem, raynum) = dynam.dpdB(:,elem-1,raynum) - ...
                gk*dot(dynam.dpdB(:,elem-1,raynum), r(:,segNum)*Vgr(segNum))/dot(gk,r(:,segNum)*Vgr(segNum));
            dynam.dxdA(:,elem,raynum) = dynam.dxdA(:,elem-1,raynum) + ...
                (Vgr(segNum)*r(:,segNum)-Vgr(segNum-1)*r(:,segNum-1))*dot(dynam.dxdA(:,elem-1,raynum),gk)/dot(Vgr(segNum-1)*r(:,segNum-1),gk);
            dynam.dxdB(:,elem,raynum) = dynam.dxdB(:,elem-1,raynum) + ...
                (Vgr(segNum)*r(:,segNum)-Vgr(segNum-1)*r(:,segNum-1))*dot(dynam.dxdB(:,elem-1,raynum),gk)/dot(Vgr(segNum-1)*r(:,segNum-1),gk);
        end
    end
    dynam.L(raynum) = sqrt( abs(det([dynam.dxdA(:,end,raynum) dynam.dxdB(:,end,raynum) Vgr(end)*r(:,end)])) / Vgr(end) / (n(2,1)/sqrt(n(1,1)^2+n(2,1)^2)) );

    %% reflection/transmission coefficients
    for RTnum = 1:noRT
        %wavetype  = raydata(1,RTnum);
        layer     = raydata(2,RTnum);
        layerNext = raydata(2,RTnum+1);
        if RTnum>1
            layerPrev = raydata(2,RTnum-1);
        else
            layerPrev = 0;
        end
        % medium 1
        AijR = Aij(:,:,layer);
        rhoR = rho(layer);
        % medium 2
        if layer == layerNext
            if layer+layerNext == 2
                AijT = Aij(:,:,2);
                rhoT = rho(2);
            else
                AijT = Aij(:,:,layer+(layer-layerPrev));
                rhoT = rho(layer+(layer-layerPrev));
            end
        else
            AijT = Aij(:,:,layerNext);
            rhoT = rho(layerNext);
        end
        
        [RT,v2,U2,P2] = RTcoef(AijR,AijT,rhoR,rhoT,n(:,RTnum),[0;0;1],U(:,RTnum));
        dynam.RTcoef(RTnum,raynum) = RT( raydata(1,RTnum+1) + abs(layerNext-layer)*3 );
    end
    
    %% calculating amplitudes
    for RTnum = 1:noRT
        dynam.P(raynum) = dynam.P(raynum) + abs(dynam.RTcoef(RTnum,raynum))*r(3,RTnum)/r(3,RTnum+1)*...
            sqrt(rho(layer)*Vgr(RTnum)/rho(layerNext)/Vgr(RTnum+1));
    end
    dynam.Amp(raynum) = dynam.P(raynum)/dynam.L(raynum)*sqrt(rho(raydata(2,1))*Vgr(1)/rho(raydata(2,end))/Vgr(end));
end
clear p3;

end    % of the function

%%