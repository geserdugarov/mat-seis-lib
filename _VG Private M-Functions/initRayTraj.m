%% *initRayTraj*
% Initialize ray trajectory

%%
% *Input:*

% xSou     - [3, 1] source coordinates
% xRec     - [3, 1] receiver coordinates 
% iRefl    - [scalar] number of the reflecting interface  
%            (*) iRefl = 0 denotes a direct wave
% zInt     - [3, :] z-coefficients of interfaces produced by 'planeEquation'
% fRefl    - [scalar] number of the reflecting fault
% zFlt     - [3, :] z-coefficients of faults produced by 'planeEquation'
% fltLayer - [scalar] layer number in which the fault reflection takes place 
% rayCode  - [3, noSeg] array resulting from 'produceRayCode'

%%
% *Output:*

% xy       - [2*(noSeg-1), 1] array containing pairs of the x- and y-coordinates of 
%            intersections of the ray trajectory with model interfaces 

%%
% *Author:* Vladimir Grechka 2012 

%%
function [xy] = initRayTraj(xSou, xRec, iRefl, zInt, fRefl, zFlt, fltLayer, rayCode)
%% Settings  
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

noSeg = size(rayCode, 2);                   % the number of ray segments
noInt = size(zInt, 2);                      % the number of model interfaces   

%% Main loop
if noSeg == 1
    xy = zeros(0, 1);                       % trajectory contains a single segment, xy is empty
else
    xy = zeros(2*(noSeg-1), 1);             % trajectory has more than one segment
    ray1 = xSou;               

    if iRefl == 0  &&  fRefl == 0
        ray2 = xRec;                        % direct ray                

    elseif iRefl ~= 0  &&  fRefl == 0       % reflection from interface iRefl
        ray2 = (xSou(1:2) + xRec(1:2))/2;   % approximation to the reflection point
        ray2(3) = dot(zInt(1:2,iRefl), ray2(1:2)) + zInt(3,iRefl);

    elseif iRefl == 0  &&  fRefl ~= 0       % reflection from fault fRefl
        ii = find(abs(zFlt(1:2,fRefl)) == max(abs(zFlt(1:2,fRefl))));
        ij = 3 - ii;                        % switch between the indexes 1 and 2 since a fault 
        ray2 = NaN(3,1);                    % is presumed to be steep
        ray2([ij,3]) = (xSou([ij,3],1) + xRec([ij,3],1))/2;     % no-interface solution
        ray2(ii) = (ray2(3) - zFlt(ij,fRefl)*ray2(ij) - zFlt(3,fRefl))/zFlt(ii,fRefl);
        
        if size(zInt, 2) ~= 0               % fault crossing model interfaces
            % Setup to move ray2(3) in increments
            zmin = min(zInt(3,:));     zmax = max(zInt(3,:));
            noInc = 20;   noPad = 5;  
            dz = (zmax - zmin)/noInc;    
            zz = zmin - noPad*dz : dz : zmax + noPad*dz;
            % Determine the depths of interfaces bounding 'fltLayer'
            if fltLayer == 1                % the fault segment is in the top layer
                zup = -Inf;
                zdn = dot(zInt(1:2,fltLayer), ray2(1:2)) + zInt(3,fltLayer);
            elseif fltLayer == noInt+1      % the fault segment is in the bottom layer
                zup = dot(zInt(1:2,fltLayer-1), ray2(1:2)) + zInt(3,fltLayer-1);
                zdn = Inf;
            else                            % the fault segment is in the middle of the model
                zup = dot(zInt(1:2,fltLayer-1), ray2(1:2)) + zInt(3,fltLayer-1);
                zdn = dot(zInt(1:2,fltLayer),   ray2(1:2)) + zInt(3,fltLayer);
            end;
            
            ray2(3) = NaN;
            for iz=1:length(zz)
                if zup < zz(iz)  &&  zz(iz) < zdn 
                    ray2(3) = zz(iz);
                    ray2(ii) = (ray2(3) - zFlt(ij,fRefl)*ray2(ij) - zFlt(3,fRefl))/zFlt(ii,fRefl);
                    break; 
                end;
            end;
            
            if isnan(ray2(3)) == 1
                fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
                display(rayCode)
                display(fRefl) 
                display(fltLayer)
                fprintf('>>> No satisfactory initial guess found \n \n');
                error('>>> STOP');
            end;
        end;
        
    else
        fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
        fprintf('>>> Erroneous combination of waveCode(3) = %g and waveCode(4) = %g \n', ...
                [iRefl, fRefl]);
          error('>>> STOP');
    end;
        
    %% Initialize trajectory assuming a straight ray 
    for iseg=1:(noSeg-1)

        if rayCode(4,iseg) == 0             % logic to switch from an interface to a fault
            iInt = rayCode(3,iseg);   fInt = 0;
            plane = zInt(:,iInt);
        else
            iInt = 0;   fInt = rayCode(4,iseg);
            plane = zFlt(:,fInt);
        end;
        
        dray = ray1 - ray2;
        den = dray(3) - dot(plane(1:2), dray(1:2));
        if abs(den) < 1.e-6 
            step = 0.5;
        else
            step = (dot(plane(1:2), ray2(1:2)) + plane(3) - ray2(3))/den;
        end;
        
        % x- and y-coordinates at which a straight ray passing through points ray1 and ray2
        % intersects plane iIint 
        xy(2*iseg-1 : 2*iseg) = ray1(1:2)*step + ray2(1:2)*(1 - step);
        if (iInt ~= 0  &&  iInt == iRefl)  ||  (iInt == 0  &&  fInt == fRefl)  
            ray1 = xRec;
        end;
    end;  
end;

end    % of the function