%% *setConstraints*
% Set constraints imposed by the stability conditions and layer velocities 

%%
% *Input:*

% Cij       - [6, 6, noLayer] array of the stiffness matrices of layers 
%   noLayer - the number of layers
% velBound  - [noLayer, 4] array of the minimum and maximum P- and S-wave velocities in layers
%             created in function 'setInitGuess'
% wellVec   - [3, 1] vector of the well direction along which the velocities are constrainted,
%             presumably by the available sonic logs
% factor    - [2, 1] vector of amplification factors for constraints imposed by the stability 
%             conditions and the layer velocities

%%
% *Output:*

% constrFun - [10*noLayer, 1] array of the constraint violations for 6 principal minors of
%             the stiffness matrices and for 4 velocity bounds defined by variable velBound
% constrJac - [10*noLayer, 21*noLayer] matrix of the derivatives of constrFun with 
%             respect to 21 stiffness components in the layers

%%
% *Author:* Vladimir Grechka 2012 - 2014

%%
function [constrFun, constrJac] = setConstraints(Cij, velBound, wellVec, factor)
%% Settings 
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

noLayer = size(Cij, 3);
F1 = zeros(6*noLayer, 1);   J1 = zeros(6*noLayer, 21*noLayer);
F2 = zeros(4*noLayer, 1);   J2 = zeros(4*noLayer, 21*noLayer);
[v2m, ~, ~, ~] = indexesVMT(1);  

%% Loop over the layers
for ilayer = 1:noLayer
    irow1 = 6*(ilayer - 1);   irow2 = 4*(ilayer - 1);   icol = 21*(ilayer - 1);
    
    %% Evaluate the stability conditions
    [vrfStabCond, detStabCond] = isStable(Cij(:,:,ilayer), 0); 
    if vrfStabCond == 0
        fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
        fprintf('>>> The elastic stability conditions are violated in layer %g \n', ilayer);

        %% Calculate the derivatives of minors with respect to Cij
        for ij = 1:21
            dCij = zeros(6,6);
            dCij(v2m(1,ij), v2m(2,ij)) = 1;
            dCij(v2m(2,ij), v2m(1,ij)) = 1;    
            
            for idet = 1:6
                if detStabCond(idet) < 0
                    F1((irow1 + idet), 1) = detStabCond(idet);
                    J1((irow1 + idet), (icol + ij)) = ...
                        dDetdX(Cij(1:idet, 1:idet, ilayer), dCij(1:idet, 1:idet));
                end;
            end;
        end;        
    end;  

    %% Evaluate the velocity constraints
    if sum(isnan(velBound(ilayer,:))) < 4
        wellVec0 = wellVec/norm(wellVec);  
        [V, U] = velPhaseU(Cij(:,:,ilayer), wellVec0);
        
        % Constrains on the S-wave velocities
        if isnan(velBound(ilayer,1)) == 0
            if V(3) < velBound(ilayer,1)
                % The slow S-wave velocity is smaller than VSmin
                F2((irow2 + 1), 1) = V(3)^2/velBound(ilayer,1)^2 - 1;
                [~, dVrow] = dVdCij(wellVec0, V(3), U(:,3));
                J2((irow2 + 1), (icol + 1) : (icol + 21)) = 2*V(3)*dVrow/velBound(ilayer,1)^2;
                fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
                fprintf('>>> The S-wave velocity constraint is violated in layer %g \n', ilayer);
                fprintf('>>> The velocity = %g is smaller than its lower bound = %g \n', ...
                        [V(3), velBound(ilayer,1)]);
            end;
        end;
        
        if isnan(velBound(ilayer,2)) == 0        
            if velBound(ilayer,2) < V(2) 
                % The fast S-wave velocity is greater than VSmax
                F2((irow2 + 2), 1) = V(2)^2/velBound(ilayer,2)^2 - 1;
                [~, dVrow] = dVdCij(wellVec0, V(2), U(:,2));
                J2((irow2 + 2), (icol + 1) : (icol + 21)) = 2*V(2)*dVrow/velBound(ilayer,2)^2;
                fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
                fprintf('>>> The S-wave velocity constraint is violated in layer %g \n', ilayer);
                fprintf('>>> The velocity = %g is greater than its upper bound = %g \n', ...
                        [V(2), velBound(ilayer,2)]);
            end;
        end;
            
        % Constrains on the P-wave velocities
        if isnan(velBound(ilayer,3)) == 0        
            if V(1) < velBound(ilayer,3)
                % The P-wave velocity is smaller than VPmin
                F2((irow2 + 3), 1) = V(1)^2/velBound(ilayer,3)^2 - 1;
                [~, dVrow] = dVdCij(wellVec0, V(1), U(:,1));
                J2((irow2 + 3), (icol + 1) : (icol + 21)) = 2*V(1)*dVrow/velBound(ilayer,3)^2;
                fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
                fprintf('>>> The P-wave velocity constraint is violated in layer %g \n', ilayer);
                fprintf('>>> The velocity = %g is smaller than its lower bound = %g \n', ...
                        [V(1), velBound(ilayer,3)]);
            end;
        end;
        
        if isnan(velBound(ilayer,4)) == 0        
            if velBound(ilayer,4) < V(1)
                % The P-wave velocity is greater than VPmax
                F2((irow2 + 4), 1) = V(1)^2/velBound(ilayer,4)^2 - 1;
                [~, dVrow] = dVdCij(wellVec0, V(1), U(:,1));
                J2((irow2 + 4), (icol + 1) : (icol + 21)) = 2*V(1)*dVrow/velBound(ilayer,4)^2;
                fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
                fprintf('>>> The P-wave velocity constraint is violated in layer %g \n', ilayer);
                fprintf('>>> The velocity = %g is greater than its upper bound = %g \n', ...
                        [V(1), velBound(ilayer,4)]);
            end;
        end;
    end;
    
end;    % of the loop over layers

%% Apply input factors to the constraints and their gradients
maxF1 = max(abs(F1));  maxJ1 = max(max(abs(J1)));
maxF2 = max(abs(F2));  maxJ2 = max(max(abs(J2)));
if maxF1 > 0
    F1 = factor(1)*F1/maxF1;   
end;
if maxJ1 > 0
    J1 = factor(1)*J1/maxJ1;
end;
if maxF2 > 0
    F2 = factor(2)*F2/maxF2;
end;
if maxJ2 > 0
    J2 = factor(2)*J2/maxJ2;
end;

constrFun = cat(1, F1, F2);
constrJac = cat(1, J1, J2);

end    % of the function