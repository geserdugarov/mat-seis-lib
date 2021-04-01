%% *getCij*
% Construct elastic parameters in the course of inversion

%% 
% *Input:*

% param        - [1, noParam] array of the stiffness-related parameters
%                noParam = length(find(unknInd == 0))  
% anisType     - [noLayer, 3][char] array indicating the layer symmetry:
%                .ISO - isotropy 
%                .VTI - transverse isotropy; the symmetry axis can be tilted
%                .ORT - orthotropy; the symmetry planes can be arbitrarily rotated 
%                .MNC - monoclinic medium with a horizontal symmetry plane; subsequent rotation
%                       is possible
%                .TRI - triclinic symmetry
% anisParam    - [noLayer, 3] string array whose rows are equal to 'Cij' or 'Ani' to indicate 
%                whether the unknowns are stiffnesses or anisotropy coefficients
% unknInd      - [noLayer, 21] index array for the elastic parameters produced by 'setUnknowns'
% flagStabCond - [scalar] flag equal to 0 or 1 to specify whether the stability conditions  
%                need to be checked

%%
% *Output:*

% CijRot       - [6, 6, noLayer] array of the stiffness coefficients in appropriately rotated
%                coordinates
% CijSym       - [6, 6, noLayer] array of the stiffness coefficients in the natural frames
% rotAngle     - [3, noLayer] array of three angles specifying the coordinate rotations 
%                used as input to 'loc2glb'
% rotMatrix    - [3, 3, noLayer] array of unitary matrices decribing the coordinate rotations
%                produced by 'loc2glb' 
% noLayer      - the number of layers
% vrfStabCond  - [1, noLayer] array of flags indicating whether the elastic stability conditions
%                are satisfied [vrfStabCond(1,ilayer) = 1] or not [vrfStabCond(1,ilayer) = -1]
% detStabCond  - [6, noLayer] array of six determinants of the principal minors of matrices
%                'CijSym'

%%
% *Author:* Vladimir Grechka 2012 2013

%%
function [CijRot, CijSym, rotAngle, rotMatrix, densities, noLayer, vrfStabCond, detStabCond] = ...
    getCij(param, anisType, anisParam, unknInd, flagStabCond)
%% Settings and defaults
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

%% Initialize Cij's
[CijRot, CijSym, rotAngle, rotMatrix, densities, noLayer] = setCij(anisType);
vrfStabCond = NaN(1, noLayer);
detStabCond = NaN(6, noLayer);

%% Compute the elastic properties
[v2m, ~, ~, ~] = indexesVMT(1);  
countParam = 0;   % counter of quantities in the 'param' array
for ilayer=1:noLayer
    C0 = zeros(6,6);

    if strcmp(anisParam(ilayer,:), 'Ani') == 1
        % Default anisotropy coefficients
        aniMNC = cij2grechka(CijSym(:,:,ilayer));
        
        if strcmp(anisType(ilayer,:), 'TRI') == 1
            fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
            fprintf(['>>> Anisotropy coefficients are undefined in layer %g ', ...
                     ' whose anisType = %s \n'], ilayer, anisType(ilayer,:));
            fprintf(['>>> Change anisType to ISO, VTI, ORT or MNC or switch to ', ...
                     'the investion in terms of Cij''s \n \n']);
              error('>>> STOP');
        end;
    end;
    
    %% Isotropic media
    if strcmp(anisType(ilayer,:), 'ISO') == 1
        indISO = [1 16];                                    % apply the v2m conversion rule in
        if strcmp(anisParam(ilayer,:), 'Cij') == 1          % function 'indexesVMT'
            for i=1:length(indISO);
                % Default to the initial model  
                    C0(v2m(1,indISO(i)), v2m(2,indISO(i))) = ...
                CijSym(v2m(1,indISO(i)), v2m(2,indISO(i)), ilayer);
                if unknInd(ilayer,i) == 0 
                    % Update the stiffness element 
                    countParam = countParam + 1;   
                    C0(v2m(1,indISO(i)), v2m(2,indISO(i))) = param(countParam); 
                end;  
            end;
            % Fill in the remaining elements
            C0 = thomsen2cij([sqrt(C0(1,1)), sqrt(C0(4,4)), 0, 0, 0]);

        elseif strcmp(anisParam(ilayer,:), 'Ani') == 1
            aniCoef = aniMNC(1:2); 
            for i=1:length(indISO);
                if unknInd(ilayer,i) == 0 
                    % Update the anisotropy coefficient 
                    countParam = countParam + 1;   
                    aniCoef(i) = param(countParam); 
                end;  
            end;
            C0 = thomsen2cij([aniCoef(1:2), 0, 0, 0]);
    
        else
            fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
            fprintf('>>> Incorrect value of anisParam = %s in layer %g \n \n', ...
                    anisParam(ilayer,:), ilayer);
              error('>>> STOP');
        end;

        if prod(unknInd(ilayer,3:21)) == 0 
            fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
            display(unknInd(ilayer,:))
            fprintf(strcat('>>> Array ''unknInd'' for isotropic layer %g has incorrect ', ...
                           ' zero values -- PAUSE \n', ilayer));    pause
        end;
    end;
    
    %% Transversely isotropic media
    if strcmp(anisType(ilayer,:), 'VTI') == 1
        indVTI = [1 3 12 16 21];                            % apply the v2m conversion rule in
        if strcmp(anisParam(ilayer,:), 'Cij') == 1          % function 'indexesVMT'
            for i=1:length(indVTI);
                % Default to the true model  
                    C0(v2m(1,indVTI(i)), v2m(2,indVTI(i))) = ...
                CijSym(v2m(1,indVTI(i)), v2m(2,indVTI(i)), ilayer);
                if unknInd(ilayer,i) == 0 
                    % Update the stiffness element 
                    countParam = countParam + 1;   
                    C0(v2m(1,indVTI(i)), v2m(2,indVTI(i))) = param(countParam); 
                end;  
            end;
            % Fill in the remaining elements
            C0(2,2) = C0(1,1);    C0(2,3) = C0(1,3);    C0(5,5) = C0(4,4);
            C0(1,2) = C0(1,1) - 2*C0(6,6);
            C0 = symMat(C0);

        elseif strcmp(anisParam(ilayer,:), 'Ani') == 1
            aniCoef = [aniMNC(1:3), aniMNC(5), aniMNC(8)]; 
            for i=1:length(indVTI);
                if unknInd(ilayer,i) == 0 
                    % Update the anisotropy coefficient 
                    countParam = countParam + 1;   
                    aniCoef(i) = param(countParam); 
                end;  
            end;
            C0 = thomsen2cij(aniCoef);

        else
            fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
            fprintf('>>> Incorrect value of anisParam = %s in layer %g \n \n', ...
                    anisParam(ilayer,:), ilayer);
              error('>>> STOP');
        end;
        
        % Rotation angles
        for i=1:2
            if unknInd(ilayer,length(indVTI)+i) == 0 
                % Update the rotation angle 
                countParam = countParam + 1;   
                rotAngle(i,ilayer) = param(countParam); 
            end;
        end;  
    
        if unknInd(ilayer,length(indVTI)+3) == 0
            fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
            fprintf('>>> Three rotation angles for VTI layer %g \n', ilayer);
            fprintf(['>>> The third angle will not be discarded because ', ...
                     'VTI stiffness tensors are invariant with respect to rotation \n']);
            fprintf('>>> around its symmetry axis -- PAUSE \n \n');   pause;
        end;
    end;
            
    %% Orthorhombic media
    if strcmp(anisType(ilayer,:), 'ORT') == 1
        indORT = [1 2 3 7 8 12 16 19 21];                   % apply the v2m conversion rule in
        if strcmp(anisParam(ilayer,:), 'Cij') == 1          % function 'indexesVMT'
            for i=1:length(indORT);
                % Default to the true model  
                    C0(v2m(1,indORT(i)), v2m(2,indORT(i))) = ...
                CijSym(v2m(1,indORT(i)), v2m(2,indORT(i)), ilayer);
                if unknInd(ilayer,i) == 0 
                    % Update the stiffness element 
                    countParam = countParam + 1;   
                    C0(v2m(1,indORT(i)), v2m(2,indORT(i))) = param(countParam); 
                end;  
            end;
            C0 = symMat(C0);

        elseif strcmp(anisParam(ilayer,:), 'Ani') == 1
            aniCoef = aniMNC(1:9); 
            for i=1:length(indORT);
                if unknInd(ilayer,i) == 0 
                    % Update the anisotropy coefficient 
                    countParam = countParam + 1;   
                    aniCoef(i) = param(countParam); 
                end;  
            end;
            C0 = tsvankin2cij(aniCoef);

        else
            fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
            fprintf('>>> Incorrect value of anisParam = %s in layer %g \n \n', ...
                    anisParam(ilayer,:), ilayer);
              error('>>> STOP');
        end;
        
        % Rotation angles
        for i=1:3
            if unknInd(ilayer,length(indORT)+i) == 0 
                % Update the rotation angle 
                countParam = countParam + 1;   
                rotAngle(i,ilayer) = param(countParam); 
            end;
        end;  
    end;
    
    %% Monoclinic media
    if strcmp(anisType(ilayer,:), 'MNC') == 1
        indMNC = [1 2 3 6 7 8 11 12 15 16 19 21];           % apply the v2m conversion rule in
        if strcmp(anisParam(ilayer,:), 'Cij') == 1          % function 'indexesVMT'
            for i=1:length(indMNC);
                % Default to the true model  
                    C0(v2m(1,indMNC(i)), v2m(2,indMNC(i))) = ...
                CijSym(v2m(1,indMNC(i)), v2m(2,indMNC(i)), ilayer);
                if unknInd(ilayer,i) == 0 
                    % Update the stiffness element 
                    countParam = countParam + 1;   
                    C0(v2m(1,indMNC(i)), v2m(2,indMNC(i))) = param(countParam); 
                end;  
            end;
            C0 = symMat(C0);

        elseif strcmp(anisParam(ilayer,:), 'Ani') == 1
            aniCoef = aniMNC(1:12); 
            for i=1:length(indMNC);
                if unknInd(ilayer,i) == 0 
                    % Update the anisotropy coefficient 
                    countParam = countParam + 1;   
                    aniCoef(i) = param(countParam); 
                end;  
            end;
            C0 = grechka2cij(aniCoef);

        else
            fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
            fprintf('>>> Incorrect value of anisParam = %s in layer %g \n \n', ...
                    anisParam(ilayer,:), ilayer);
              error('>>> STOP');
        end;
        
        % Rotation angles
        for i=1:3
            if unknInd(ilayer,length(indMNC)+i) == 0 
                % Update the rotation angle 
                countParam = countParam + 1;   
                rotAngle(i,ilayer) = param(countParam); 
            end;
        end;  
    end;
        
    %% Triclinic media
    if strcmp(anisType(ilayer,:), 'TRI') == 1
        indTRI = 1:21;  
        for i=1:length(indTRI);
            % Default to the true model  
                C0(v2m(1,indTRI(i)), v2m(2,indTRI(i))) = ...
            CijSym(v2m(1,indTRI(i)), v2m(2,indTRI(i)), ilayer);
            if unknInd(ilayer,i) == 0 
                % Update the stiffness element 
                countParam = countParam + 1;   
                C0(v2m(1,indTRI(i)), v2m(2,indTRI(i))) = param(countParam); 
            end;  
        end;
        C0 = symMat(C0);
    end;
            
    %% Check the stability conditions
    if flagStabCond == 1
        [vrfStabCond(1,ilayer), detStabCond(:,ilayer)] = isStable(C0, 0); 
        if vrfStabCond(1,ilayer) < 0;
            fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
            fprintf('>>> The elastic stability conditions are violated in layer %g for tensor', ...
                ilayer);
            display(C0);
        end;  
    end;

    %% Save stiffness tensor in the natural coordinate frame
    CijSym(:,:,ilayer) = C0; 
    
    %% Apply rotation  
    if isempty(find(rotAngle, 1)) == 0
        % Rotate the stiffness matrix
        rotMatrix(:,:,ilayer) = loc2glb(rotAngle(:,ilayer));   
        [CijRot(:,:,ilayer), ~] = bond(C0, rotMatrix(:,:,ilayer));
    else
        CijRot(:,:,ilayer) = C0;
    end;
            
end;   % of loop over layers

end    % of the function
