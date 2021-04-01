%% *reduceDerivMat*
% Tailor a full matrix of the Frechet derivatives of traveltimes to a given inverse problem 

%%
% *Input:*

% inpDT            - Frechet-derivative matrix of the traveltimes 
%                    For 'inpDT' produced by 'art3D',
%                    size(inpDT,1) = noTime  = noSou*noRec*noWave,
%                    size(inpDT,2) = noParam = 21*(noInt-1) + 3*noSou + ...
%                                                             3*noRec + 3*noInt + 3*noFlt + noSou 
% dim              - dimension of the event-location problem:
%                    . dim = 2 for 2D and 
%                    . dim = 3 for 3D and 
%
% model.           - model structure whose following fields are used:
%     Cij.symmetry - [6, 6, noLayer] array of the stiffness matrices in crystallographic 
%                    coordinates
%     rotation.    - nested structure with two fields:
%         angles   - [3, noLayer] array of three angles [azim, tilt, azx1] described in 'loc2glb'
%                    that specify the coordinate rotations in layers
%         matrix   - [3, 3, noLayer] array of unitary matrices decribing the coordinate rotations
%                    typically resulting from 
%                    'rotation.matrix(:,:,ilayer) = loc2glb(rotation.angles(:,ilayer))'
%     anisType     - [noLayer, 3][char] array indicating the layer symmetry:
%                    . ISO - isotropy 
%                    . VTI - transverse isotropy; the symmetry axis can be tilted
%                    . ORT - orthotropy; the symmetry planes can be arbitrarily rotated 
%                    . MNC - monoclinic medium with a horizontal symmetry plane; subsequent rotation
%                            is possible
%                    . TRI - triclinic symmetry
%     anisParam    - [noLayer, 3] character array whose rows are equal to 'Cij' or 'Ani' to indicate 
%                    whether the unknowns are stiffnesses or anisotropy coefficients
%     azimSou      - [1, noSou] array of the azimuths of sources from the receiver string
%                    (*) If dim = 3, azimSou is irrelevant
%
% unknowns.        - structure of unknowns whose following fields are used:
%     CijInd       - [noLayer, 21] index array for the elastic parameters
%     interfaceInd - [5, noInt] index array for the interfaces
%     faultInd     - [5, noFlt] index array for the faults
%     perfInd      - [1, noSou] index array for the perforations
%     tauInd       - [1, noSou] index array for the origin times
%                  
% (*) Array sizes:
%          noLayer - the number of layers
%          noSou   - the number of sources
%          noRec   - the number of receivers
%          noInt   - the number of interfaces
%          noFlt   - the number of faults

%%
% *Output:*

% outDT            - [noTime, noUnkn] Frechet-derivative matrix containing columns corresponding 
%                    to the unknown parameters defined in 'setUnknowns'
%                    noUnkn = length(find(unknowns.CijInd == 0)) + ...
%                             length(find(unknowns.interfaceInd == 0)) + ...
%                             length(find(unknowns.faultInd == 0)) + ...
%                         dim*length(find(unknowns.perfInd == 0)) + ...
%                             length(find(unknowns.tauInd == 0))

%%
% *Author:* Vladimir Grechka 2012 - 2014

%%
function [outDT] = reduceDerivMat(inpDT, dim, model, unknowns)
%% Settings 
[~, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

% Unpack some of the 'model' structure
noSou = model.noSou;   noRec = model.noRec;   noInt = model.noInt;   noFlt = model.noFlt;    
anisType = model.anisType;           anisParam = model.anisParam;

noLayer = size(anisType, 1);
[noTime, noParam] = size(inpDT);
noUnkn = length(find(unknowns.CijInd == 0)) + length(find(unknowns.interfaceInd == 0)) + ...
         length(find(unknowns.faultInd == 0)) + dim*length(find(unknowns.perfInd == 0)) + ...
         length(find(unknowns.tauInd == 0));
outDT = NaN(noTime, noUnkn);

if ~(dim == 2  ||  dim == 3)   
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Incorrect dimension of the even-location problem dim = %g \n', dim);
    fprintf('>>> The correct values are dim = 2 or dim = 3 \n \n');  
      error('>>> STOP');
end;

%% Derivatives with respect to elastic parameters of layers
inpCount = 0;   outCount = 0;
for ilayer = 1:noLayer
    
    % Set the maximum numbers of unknown stiffness or anisotropy coefficients
    if strcmp(anisType(ilayer,:), 'ISO') == 1;   noUnknCij = 2;   end;
    if strcmp(anisType(ilayer,:), 'VTI') == 1;   noUnknCij = 5;   end;
    if strcmp(anisType(ilayer,:), 'ORT') == 1;   noUnknCij = 9;   end;
    if strcmp(anisType(ilayer,:), 'MNC') == 1;   noUnknCij = 12;  end;
    if strcmp(anisType(ilayer,:), 'TRI') == 1;   noUnknCij = 21;  end;
    
    %----------------------------------------------------------------------------------------------
    % Derivatives with respect to the stiffness or anisotropy coefficients    
    clear indUnkn inpDTSym    
    indUnkn = find(unknowns.CijInd(ilayer,1:noUnknCij) == 0);
    count = length(indUnkn);

    if count > 0    
        % Rotate derivatives d/dCij to the natural coordinate frame in which the unknowns are 
        % specified 
        inpDTSym = dFdCijRot(inpDT(:, (inpCount + 1) : (inpCount + 21)), ...
                             model.rotation.matrix(:,:,ilayer)');
        
        if strcmp(anisParam(ilayer,:), 'Cij') == 1
            % Derivatives d/dCij for symmetries higher than triclinic
            [outDT(:, (outCount + 1) : (outCount + count))] = ...
                dFdCijSym(inpDTSym, anisType(ilayer,:), indUnkn);
        
        elseif strcmp(anisParam(ilayer,:), 'Ani') == 1
            % Frechet derivatives with respect to the anisotropy parameters
            [outDT(:, (outCount + 1) : (outCount + count)), ~] = ...
                dFdAniCoef(inpDTSym, model.Cij.symmetry(:,:,ilayer), anisType(ilayer,:), indUnkn);
        else
            fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
            fprintf('>>> Variable ''anisParam'' in layer %g is %s \n', ...
                    ilayer, anisParam(ilayer,:));
            fprintf('>>> Its correct values are ''Ani'' or ''Cij'' \n \n');
              error('>>> STOP');
        end;        
        outCount = outCount + count;
    end;

    %----------------------------------------------------------------------------------------------
    % Derivatives with respect to the rotation angles    
    if strcmp(anisType(ilayer,:), 'TRI') == 0
        clear indUnknRot     
        indUnknRot = find(unknowns.CijInd(ilayer, (noUnknCij+1) : (noUnknCij+3)) == 0);
        count = length(indUnknRot);
    
        if count > 0         
            [outDT(:, (outCount + 1) : (outCount + count))] = ...
                dFdRotAngles(inpDT(:, (inpCount + 1) : (inpCount + 21)), ...
                             model.Cij.symmetry(:,:,ilayer), ...
                             model.rotation.angles(:,ilayer), indUnknRot);                     
            outCount = outCount + count;
        end;
    end;
    
    inpCount = inpCount + 21;
            
end;  % of loop with respect to the stiffness derivatives

%% Derivatives with respect to the source coordinates
if inpCount < noParam 
    for isou = 1:noSou
        if unknowns.perfInd(isou) == 0
            if dim == 2
                outDT(:, outCount + 1) = inpDT(:, inpCount + 1)*cos(model.azimSou(isou)) + ... 
                                         inpDT(:, inpCount + 2)*sin(model.azimSou(isou)); 
                outDT(:, outCount + 2) = inpDT(:, inpCount + 3);
            else
                for i=1:3
                    outDT(:, outCount + i) = inpDT(:, inpCount + i);
                end;
            end;
            outCount = outCount + dim;
        end;
        inpCount = inpCount + 3;
    end;  % of loop with respect to the source-location derivatives
    
    if inpCount < noParam 
        % Account for the derivatives with respect to the receiver coordinates
        inpCount = inpCount + 3*noRec;  
        
        %% Derivatives with respect to the interface parameters
        [irow1, icol1] = find(unknowns.interfaceInd == 0);
        for ii1 = 1:length(irow1)
            outDT(:, outCount + ii1) = ...
            inpDT(:, inpCount + 3*(icol1(ii1)-1) + (irow1(ii1)-2));
        end;
        inpCount = inpCount + 3*noInt;  
        outCount = outCount + length(irow1); 
        
        %% Derivatives with respect to the fault parameters
        [irow2, icol2] = find(unknowns.faultInd == 0);
        for ii2 = 1:length(irow2)
            outDT(:, outCount + ii2) = ...
            inpDT(:, inpCount + 3*(icol2(ii2)-1) + (irow2(ii2)-2));
        end;
        inpCount = inpCount + 3*noFlt;  
        outCount = outCount + length(irow2); 

        if inpCount < noParam
            %% Derivatives with respect to the origin times
            count = 0;
            for isou = 1:noSou
                if unknowns.tauInd(isou) == 0
                    count = count + 1;
                    outDT(:, outCount + count) = inpDT(:, inpCount + isou);
                end;
            end;   

        end;
    end;
end;

end    % of the function