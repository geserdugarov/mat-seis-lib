%% *rayAttributes*
% Collect attributes along a ray trajectory

%%
% *Input:*

% xy         - [2*(noSeg-1), 1] array of the pairs of x- and y-coordinates of intersections of a 
%              given ray with the model interfaces and faults
% rayType    - [noSeg, 1] array of ray types (1 for P, 2 for S1 or SV, and 3 for S2 or SH wave) 
%              along each ray segment
%   noSeg    - the number of ray segments
% xSou       - [3, 1] source coordinates
% xRec       - [3, 1] receiver coordinates 
% zInt       - [3, noSeg-1] array of the coefficients of interfaces given by 
%              z = x*zInt(1) + y*zInt(2) + zInt(3)
%              (*) The interface depth at the i-th ray intresection is 
%                  z = dot(zInt(1:2,i), xy(2*i-1 : 2*i, 1)) + zInt(3,i)
% Cij        - [6, 6, noSeg] stiffness matrices of layers arranged along the ray trajectory
% flagWAA    - flag indicating whether calculations should be performed for weak
%              (flagWAA = 'WA') or strong (flagWAA = 'SA') anisotropy  

%%
% *Output:*

% xyz        - [noSeg+1, 3] ray-trajectory array including the source and the receiver 
% tOut       - [1, noSeg] array of times along ray segments 
% nOut       - [3, noSeg] array of the wave normals 
% rOut       - [3, noSeg] array of the unit ray directions 
% pOut       - [3, noSeg] array of the slowness vectors
% UOut       - [3, noSeg] array of the polarization vectors 
% VOut       - [1, noSeg] array of the (scalar) phase velocities 
% dVOut      - [3, noSeg] array of derivatives of the phase velocities with respect to the ray
%              direction  
% gOut       - [1, noSeg] array of the (scalar) group velocities 
% dgOut      - [2, noSeg] array of derivatives of the group velocities computed in 'velRay' 
% KpOut      - [1, noSeg] array of the Gaussian curvatures of the slowness surfaces
% rayTypeOut - [noSeg, 1] array 'rayType' corrected in such a way that that the S-wave 
%              polarization vector at ray segment 'iseg' is as close as possible to the  
%              polarization vector at segment 'iseg - 1'
%              (*) This adjustment is needed to prevent the nearly zero reflection and transmission 
%                  coefficients of the S-waves caused by the rapid polarization change 
%                  (e.g., from SV to SH or from to SH) in the adjacent layers 
% flagOut    - convergence flag produced by function velRay:
%              flagOut = 1 - solution is found
%              flagOut = 0 - no solution is found; all other output parameters are set to zero 

%%
% *Author:* Vladimir Grechka 2012 - 2014

%%
function [xyz, tOut, nOut, rOut, pOut, UOut, VOut, dVOut, gOut, dgOut, KpOut, rayTypeOut, flagOut] = ...
          rayAttributes(xy, rayType, xSou, xRec, zInt, Cij, flagWAA)
%% Settings 
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

noSeg = size(rayType, 1);
noIter = 50;    tolParam = 1.e-10;      % parameters governing iterations in function velRay

tOut  = zeros(1, noSeg);  VOut  = tOut;  gOut = tOut;  KpOut = tOut; 
dgOut = zeros(2, noSeg);    
nOut  = zeros(3, noSeg);  rOut  = nOut;  pOut = nOut;  UOut  = nOut;  dVOut = nOut;
rayTypeOut = rayType;     rayTypeFlip = [1, 3, 2];    tolFlip = 1.e-6;
flagOut = 1;

%% Coordinates of the ray trajectory
xyz = [xSou'; zeros(noSeg - 1, 3); xRec'];
for iseg = 2:noSeg
    xyz(iseg,1:2) = xy(2*iseg-3 : 2*iseg-2, 1);
    xyz(iseg,  3) = dot(xyz(iseg, 1:2), zInt(1:2, iseg-1)) + zInt(3, iseg-1);
end

%% Compute the attributes
for iseg = 1:noSeg
    raySeg = (xyz(iseg+1,:) - xyz(iseg,:))';    
    rayLength = norm(raySeg);   
    rOut(:,iseg) = raySeg/rayLength;
    [Vtmp, Utmp] = velPhaseU(Cij(:,:,iseg), rOut(:,iseg));
    
    % old: Change 'rayTypeOut' for the S-waves in isotropic layer if necessary
    % Change polarization sequence for the S-waves in isotropic layer if necessary
    if iseg > 1  &&  rayTypeOut(iseg) > 1  &&  isISO(Cij(:,:,iseg), []) == 1
        if  abs(dot( Usave(:,rayTypeOut(iseg-1)), Utmp(:,rayTypeFlip(rayTypeOut(iseg))) )) > ...
            abs(dot( Usave(:,rayTypeOut(iseg-1)), Utmp(:,rayTypeOut(iseg)) )) + tolFlip
            % Flip the shear-wave type
            Vtmp = Vtmp(rayTypeFlip);
            Utmp = Utmp(:,rayTypeFlip);
            % old: rayTypeOut(iseg) = rayTypeFlip(rayTypeOut(iseg));
        end
    end
    
    if strcmp(flagWAA, 'WA') == 1
        nOut(:,iseg) = rOut(:,iseg);
        VOut(1,iseg) = Vtmp(rayTypeOut(iseg));
        
        % Make sure that the phase velocity is real-valued
        if isreal(VOut(1,iseg)) == 0
            flagOut = 0;    
            return;    
        end 
        gOut(1,iseg) = VOut(iseg);
        UOut(:,iseg) = Utmp(:,rayTypeOut(iseg));
        
        % Compute derivatives of the phase velocity 
        g = velGroup(Cij(:,:,iseg), nOut(:,iseg), Vtmp, Utmp, rayTypeOut(iseg));
        dVOut(:,iseg) = g(:,rayTypeOut(iseg)) - nOut(:,iseg)*VOut(1,iseg);
        
        % WAA representation of the slowness vector
        pOut(:,iseg) = dTdR((xyz(iseg+1,:) - xyz(iseg,:))', VOut(iseg), dVOut(:,iseg), [], flagWAA);

        % Compute the Gaussian curvature of the slowness surface
        KpOut(1,iseg) = VOut(1,iseg)^2;
        % instead of the exact computatation
        % [n1, V1, U1, g1, dg1, KpOut(iseg), misfitRay1, flagOut1] = ...
        %     velRay1989(Cij(:,:,iseg), rOut(:,iseg), nOut(:,iseg), rayType(iseg), tolParam, 1);
    else
        % flagWAA = 'SA'
        [nOut(:,iseg), VOut(1,iseg), UOut(:, iseg), g_tmp, dgOut(:,iseg), ...
            KpOut(1,iseg), flagOut] = ...
            velRay(Cij(:,:,iseg), rOut(:,iseg), rayTypeOut(iseg), tolParam, noIter, 0);
        pOut(:,iseg) = nOut(:,iseg)/VOut(1,iseg);
        gOut(1,iseg) = norm(g_tmp);                          
        if flagOut == 0;    return;    end
    end
    

    [~, Utmp] = velPhaseU(Cij(:,:,iseg), nOut(:,iseg));
    
    % Change polarization directions for the waves if necessary
    % changing directions of two vectors for saving right-handed coordinate system
    if rayTypeOut(iseg) < 3
        secondwavetype = rayTypeOut(iseg)+1;
    else
        secondwavetype = rayTypeOut(iseg)-1;
    end
    if iseg == 1
        % in the direction of the ray
        if abs(dot( rOut(:,iseg), Utmp(:,rayTypeOut(iseg)) )) > 1.0e-6
            if  dot( rOut(:,iseg), Utmp(:,rayTypeOut(iseg)) ) < 0
                Utmp(:,rayTypeOut(iseg)) = -Utmp(:,rayTypeOut(iseg));
                Utmp(:,secondwavetype)   = -Utmp(:,secondwavetype);
            end
        else
            % for SH we can analyse qSV and change both directions
            if abs(dot( rOut(:,iseg), Utmp(:,secondwavetype) )) > 1.0e-6
                if  dot( rOut(:,iseg), Utmp(:,secondwavetype) ) < 0
                    Utmp(:,rayTypeOut(iseg)) = -Utmp(:,rayTypeOut(iseg));
                    Utmp(:,secondwavetype)   = -Utmp(:,secondwavetype);
                end
            end
        end
    else
        if sign(rOut(3,iseg-1)) == sign(rOut(3,iseg)) % before reflection
            if rayTypeOut(iseg-1) == rayTypeOut(iseg)
                if  dot( Usave(:,rayTypeOut(iseg-1)), Utmp(:,rayTypeOut(iseg)) ) < 0
                    Utmp(:,rayTypeOut(iseg)) = -Utmp(:,rayTypeOut(iseg));
                    Utmp(:,secondwavetype)   = -Utmp(:,secondwavetype);
                end
            else
                fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
                error('>>> No correct checking of polarization for transmission with changing of wave type.');
            end
        else % for reflection
            if rayTypeOut(iseg) == 1 % P-wave
                if  dot( rOut(:,iseg), Utmp(:,rayTypeOut(iseg)) ) < 0
                    Utmp(:,rayTypeOut(iseg)) = -Utmp(:,rayTypeOut(iseg));
                    Utmp(:,secondwavetype)   = -Utmp(:,secondwavetype);
                end
            else
                % search SH waves before and after reflection
                temp = [abs(dot( rOut(:,iseg), Utmp(:,1) )) abs(dot( rOut(:,iseg), Utmp(:,2) )) ...
                        abs(dot( rOut(:,iseg), Utmp(:,3) ))];
                [~,pos2] = min(temp);
                temp = [abs(dot( rOut(:,iseg-1), Usave(:,1) )) abs(dot( rOut(:,iseg-1), Usave(:,2) )) ...
                        abs(dot( rOut(:,iseg-1), Usave(:,3) ))];
                [~,pos1] = min(temp);
                % SH wave should not change direction
                if  dot( Usave(:,pos1), Utmp(:,pos2) ) < 0
                    Utmp(:,rayTypeOut(iseg)) = -Utmp(:,rayTypeOut(iseg));
                    Utmp(:,secondwavetype)   = -Utmp(:,secondwavetype);
                end
            end
        end
    end
    
    UOut(:, iseg) = Utmp(:,rayTypeOut(iseg));
    
    Usave = Utmp;   % save the polarization matrix for use at the nect ray segment
    tOut(1,iseg) = rayLength/gOut(1,iseg);
end

end    % of the function
