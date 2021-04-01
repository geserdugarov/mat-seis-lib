%% *produceRayCode*
% Make the ray code for a direct or reflected ray that connects a given source-receiver pair

%%
% *Input:*

% xSou      - [3, 1] source coordinates
% xRec      - [3, 1] receiver coordinates 
% zInt      - [3, :] array of the z-coefficients of interfaces produced by 'planeEquation'
% zFlt      - [3, :] array of the z-coefficients of faults produced by 'planeEquation'
% nrmInt    - [3, :] array of the components of normals to interfaces produced by 'planeEquation'
% nrmFlt    - [3, :] array of the components of normals to faults produced by 'planeEquation'
% Cij       - [6, 6, :] array of the stiffness matrices of layers produced by 'setCij'
% waveCode  - input wave code produced by 'setWaveCode'
%             . waveCode(1, i) = 1, 2, or 3 - for P-, S1-, or S2-wave leaving a source 
%             . waveCode(2, i) = 0, 1, 2, or 3 - for reflected P-, S1-, or S2-wave
%             . waveCode(2, i) = 0 means no reflection, hence, traveltime of a direct wave 
%                                  will be computed                            
%             . waveCode(3, i) = 1, ..., noInt - the number of reflecting interface 
%             . waveCode(4, i) = 1, ..., noFault - the number of reflecting fault 
%             (*) If waveCode(2, i) = 0, the values of waveCode([3,4], i) are irrelevant   
% fltLayer  - the number of layer containing the examined fault segment 

%%
% *Output:*

% rayCode   - [4, noSeg] array containing the rows:
%             . rayCode(1, :) - ray type = 1, 2, or 3 for the P-, S1-, or S2-wave along each 
%                               ray segment
%             . rayCode(2, :) - the sequence of layer numbers sequentially crossed by the ray
%             . rayCode(3, :) - the sequence of interface numbers
%             . rayCode(3, :) = 0 indicates reflection from a fault
%             . rayCode(4, :) - the sequence of numbers of reflecting faults 
%             . rayCode(4, :) = 0 indicates the absence of reflection from a fault at a given
%                                 ray segment
%     noSeg - the number of ray segments
% sou3      - the vertical source coordinate
%             (*) If necessary, it is perturbed to move the source away from the model interface
% rec3      - the vertical receiver coordinate
%             (*) If necessary, it is perturbed to move the receiver away from the model interface
% zIntRay   - [3, noSeg-1] array of the z-coefficients of interfaces and faults intersected
%             by the ray trajectory
% nrmIntRay - [3, noSeg-1] array of the components of normals to interfaces and faults intersected
%             by the ray trajectory
% CijRay    - [6, 6, noSeg] array of the stiffness matrices of layers along the trajectory

%%
% *Author:* Vladimir Grechka 2012 - 2014

%%
function [rayCode, sou3, rec3, zIntRay, nrmIntRay, CijRay] = ...
        produceRayCode(xSou, xRec, zInt, zFlt, nrmInt, nrmFlt, Cij, waveCode, fltLayer)
%% Settings 
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

depthPert = 1.e-6;                              % perturbation of the source and receiver depths 
sou3 = xSou(3);                                 % intended to move them away from an interface
rec3 = xRec(3); 
souLayer = [];   recLayer = [];
noInt = size(zInt,2);   noFlt = size(zFlt,2);

%% Determine the numbers of layers containing the source and receiver 

% Perturb the source and receiver depths, if necessary 
if noInt > 0
    zIntSou = dot(zInt, repmat([xSou(1:2,1); 1], 1, noInt), 1);
    zIntRec = dot(zInt, repmat([xRec(1:2,1); 1], 1, noInt), 1);
end; 
 
if noInt == 0
    souLayer = 1;                                               % homogeneous space
    recLayer = 1;
else
    for i=1:noInt
        if abs(zIntSou(i) - xSou(3)) < depthPert 
            sou3 = xSou(3)*(1 + depthPert);                     % perturb the source depth
            xSou(3) = sou3;
        end;
        % Determine the layer number in which the source is located
        [~, iis] = sort([xSou(3), zIntSou], 'ascend'); 
        souLayer = find(iis == 1); 

        if abs(zIntRec(i) - xRec(3)) < depthPert 
            rec3 = xRec(3)*(1 + depthPert);                     % perturb the receiver depth
            xRec(3) = rec3;
        end;
        % Determine the layer number in which the receiver is located
        [~, iir] = sort([xRec(3), zIntRec], 'ascend'); 
        recLayer = find(iir == 1);
    end;
end;

if isempty(souLayer) == 1
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Unassigned layer number for the source with coordinates [%g, %g, %g] \n', xSou');
    fprintf('>>> Please check the input \n \n');
      error('>>> STOP');
end;        
        
if isempty(recLayer) == 1
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Unassigned layer number for the receiver with coordinates [%g, %g, %g] \n', xRec');
    fprintf('>>> Please check the input \n \n');
      error('>>> STOP');
end;        

%% Check the number of input interfaces 
if waveCode(3) > noInt
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> The ''interface'' array is too small: \n');
    fprintf('        while only %g interfaces are specified, \n', noInt);
    fprintf('        the ray trajectory [%g %g %g %g] refers to the interface number %g \n', ...
            [waveCode', waveCode(3)]);
    fprintf('>>> Please increase the size of ''interface'' array or modify ''waveCode'' \n \n');
      error('>>> STOP');
end;    

%% Check the number of input faults  
if waveCode(4) > noFlt
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> The ''fault'' array is too small: \n');
    fprintf('        while only %g faults are specified, \n', noFlt);
    fprintf('        the ray trajectory [%g %g %g %g] refers to the fault number %g \n', ...
            [waveCode', waveCode(4)]);
    fprintf('>>> Please increase the size of ''fault'' array or modify ''waveCode'' \n \n');
      error('>>> STOP');
end;    

%% Determine the type of a given ray trajectory

% Find out whether a given trajectory corresponds to a direct ray (waveIndex = 1), a reflection 
% from an interface (waveIndex = 2), or a reflection from a fault (waveIndex = 3) 
if waveCode(2) == 0
    waveIndex = 1;                                              % direct ray
elseif waveCode(2) ~= 0  &&  waveCode(3) ~= 0  &&  waveCode(4) == 0
    if (dot(zInt(:,waveCode(3)), [xSou(1:2,1); 1]) - xSou(3))* ... 
       (dot(zInt(:,waveCode(3)), [xRec(1:2,1); 1]) - xRec(3)) > 0
       waveIndex = 2;                                           % reflection from an interface
    else
       waveIndex = 1;                                           % direct ray    
    end;
elseif waveCode(2) ~= 0  &&  waveCode(3) == 0  &&  waveCode(4) ~= 0  
    if (dot(zFlt(:,waveCode(4)), [xSou(1:2,1); 1]) - xSou(3))* ... 
       (dot(zFlt(:,waveCode(4)), [xRec(1:2,1); 1]) - xRec(3)) > 0
       waveIndex = 3;                                           % reflection from a fault
    else
       waveIndex = 1;                                           % direct ray    
    end;
elseif waveCode(2) ~= 0  &&  waveCode(3) ~= 0  &&  waveCode(4) ~= 0  
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Inconsistent values of ''waveCode'' = [%g %g %g %g] \n', waveCode); 
    fprintf('>>> Either ''waveCode(3)'' or ''waveCode(4)'' should be zero \n'); 
    fprintf('>>> Please correct ''waveCode'' \n \n');
      error('>>> STOP');
else
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf(['>>> The values of ''waveCode'' = [%g %g %g %g] imply a reflection from ', ...
             ' a nonexistent interface or fault \n'], waveCode); 
    fprintf('>>> Please correct ''waveCode'' \n \n');
      error('>>> STOP');
end;

%% Fill array 'rayCode'
switch waveIndex    
    case 1
        % Direct ray
        rayCodeLength = abs(recLayer - souLayer) + 1; 
        rayCode = NaN(4, rayCodeLength);
        rayCode(1,:) = waveCode(1); 
        if souLayer <= recLayer
            rayCode(2,:) = [souLayer : 1 : recLayer]; 
            rayCode(3,:) = rayCode(2,:);  
        else
            rayCode(2,:) = [souLayer : -1 : recLayer];     
            rayCode(3,:) = rayCode(2,:) - 1;  
        end;
        rayCode(4,:) = zeros(1, rayCodeLength);

    case 2
        % Reflection from an interface
        if souLayer <= waveCode(3)
            % The source is shallower than the reflector 
            % (the receiver should be shallower than the reflector too)
            rayCodeLengthInc = abs(waveCode(3) - souLayer) + 1; 
            rayCodeLengthRef = abs(recLayer - waveCode(3)) + 1; 
            rayCode = NaN(4, rayCodeLengthInc + rayCodeLengthRef);
            rayCode(2,:) = [souLayer : 1 : waveCode(3), waveCode(3) : -1 : recLayer]; 
            rayCode(3,:) = [souLayer : 1 : waveCode(3), waveCode(3)-1 : -1 : recLayer-1];
        else
            % The source is deeper than the reflector 
            % (the receiver should be deeper than the reflector too)
            rayCodeLengthInc = abs(waveCode(3) - souLayer); 
            rayCodeLengthRef = abs(recLayer - waveCode(3)); 
            rayCode = NaN*ones(4, rayCodeLengthInc + rayCodeLengthRef);
            rayCode(2,:) = [souLayer : -1 : waveCode(3)+1, waveCode(3)+1 : 1 : recLayer]; 
            rayCode(3,:) = [souLayer-1 : -1 : waveCode(3), waveCode(3)+1 : 1 : recLayer];
        end
        rayCode(1,:) = [waveCode(1)*ones(1, rayCodeLengthInc), ...
                        waveCode(2)*ones(1, rayCodeLengthRef)]; 
        rayCode(4,:) = zeros(1, (rayCodeLengthInc + rayCodeLengthRef));
                                            
    case 3
        % Reflection from a fault
        rayCodeLengthInc = abs(fltLayer - souLayer) + 1; 
        rayCodeLengthRef = abs(recLayer - fltLayer) + 1;
        rayCode = NaN(4, rayCodeLengthInc + rayCodeLengthRef);
                
        if souLayer <= fltLayer
            rayCode(2, 1:rayCodeLengthInc) = [souLayer : 1 : fltLayer]; 
            rayCode(3, 1:rayCodeLengthInc) = rayCode(2, 1:rayCodeLengthInc);  
        else
            rayCode(2, 1:rayCodeLengthInc) = [souLayer : -1 : fltLayer];     
            rayCode(3, 1:rayCodeLengthInc) = rayCode(2, 1:rayCodeLengthInc) - 1;  
        end;

        if fltLayer <= recLayer  
            rayCode(2, rayCodeLengthInc+1:end) = [fltLayer  : 1 : recLayer]; 
            rayCode(3, rayCodeLengthInc+1:end) = rayCode(2, rayCodeLengthInc+1:end);  
        else
            rayCode(2, rayCodeLengthInc+1:end) = [fltLayer  : -1 : recLayer];     
            rayCode(3, rayCodeLengthInc+1:end) = rayCode(2, rayCodeLengthInc+1:end) - 1;  
        end;
        
        rayCode(1,:) = [waveCode(1)*ones(1, rayCodeLengthInc), ...
                        waveCode(2)*ones(1, rayCodeLengthRef)]; 
        rayCode(3,rayCodeLengthInc) = 0;    
        rayCode(4,:) = zeros(1, rayCodeLengthInc + rayCodeLengthRef);
        rayCode(4,rayCodeLengthInc) = waveCode(4);
otherwise
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Error in assigning ''rayCode'' \n'); 
    display(rayCode);
      error('>>> STOP');
end;

%% Check whether the size of array 'Cij' is correct   
if max(rayCode(2,:)) > size(Cij,3)
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> The ''Cij'' array is too small: \n');
    fprintf('    ''Cij'' matrix(es) is(are) defined for %g layer(s), \n', size(Cij,3));
    fprintf('    while the ray [%g %g %g %g] intersects layer number %g \n', ...
            [waveCode', max(rayCode(2,:))]);
    fprintf('>>> Add rows to the ''anisType'' array or modify ''waveCode'' \n \n');
      error('>>> STOP');
end;    

%% Create the interface and stiffness arrays along the trajectory
noSeg = size(rayCode, 2);
zIntRay = zeros(3, noSeg-1);  nrmIntRay = zeros(3, noSeg-1);
for iseg=1:(noSeg-1)
    if rayCode(4,iseg) ~= 0
        zIntRay(:,iseg)   =   zFlt(:,rayCode(4,iseg));        % reflection from a fault
        nrmIntRay(:,iseg) = nrmFlt(:,rayCode(4,iseg));        % normal to the fault
    else 
        zIntRay(:,iseg)   =   zInt(:,rayCode(3,iseg));        % intersection with or reflection 
                                                              % from an interface
        nrmIntRay(:,iseg) = nrmInt(:,rayCode(3,iseg));        % normal to the interface
    end; 
end;

CijRay = zeros(6, 6, noSeg);
for iseg=1:noSeg
    CijRay(:,:,iseg) = Cij(:,:,rayCode(2,iseg));
end;

end    % of the function