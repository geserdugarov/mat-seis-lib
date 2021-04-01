%% *generateCaptions*
% Generate captions for various plots and determine the units of columns of the matrix of
% Frechet derivatives of travelimes

%%
% *Input:*

% dim            - [scalar] dimension of the source-location problem 
%                  dim = 2 when events are located in the vertical planes specified by the azimuths
%                  of the P-wave polarization vectors and 
%                  dim = 3 for full 3D event-location procedure  
%
% model.         - structure created by function 'setModel' whose following fields are used:
%     noSou      - the number of sources
%     anisType   - [noLayer, 3][char] array indicating the layer symmetry:
%                  .ISO - isotropy 
%                  .VTI - transverse isotropy; the symmetry axis can be tilted
%                  .ORT - orthotropy; the symmetry planes can be arbitrarily rotated 
%                  .MNC - monoclinic medium with a horizontal symmetry plane; subsequent rotation
%                         is possible
%                  .TRI - triclinic symmetry
%     anisParam  - [noLayer, 3] character array whose rows are equal to 'Cij' or 'Ani' to indicate 
%                  whether the unknowns are stiffnesses or anisotropy coefficients
%        noLayer - the number of layers
%
% unknowns.      - structure created by function 'setUnknowns' whose following fields are used:
%   CijInd       - [noLayer, 21] index array for the elastic parameters
%   interfaceInd - [5, noInt] index array for the interfaces
%          noInt - the number of interfaces
%   faultInd     - [5, noFlt] index array for the faults
%          noFlt - the number of faults
%   perfInd      - [1, noSou] index array for the perforations
%   tauInd       - [1, noSou] index array for the origin times
%          noSou - the number of sources

%%
% *Output:*

% caption        - {[':', noUnkn]} cell array containing strings that label the unknowns  
%                  noUnkn = length(find(unknowns.CijInd == 0)) + ...
%                           length(find(unknowns.interfaceInd == 0)) + ...
%                           length(find(unknowns.faultInd == 0)) + ...
%                       dim*length(find(unknowns.perfInd == 0)) + ...
%                           length(find(unknowns.tuaInd == 0))
% units          - [3, noUnkn] matrix whose rows contain exponentials of time, distance, and angle, 
%                  respectively; for example, the derivative (d time)/(d Vp0) has units [T^2/L] and  
%                  produces the column [2, -1, 0]'; likewise, the units [T/A] of the derivative 
%                 (d time)/(d tilt) yield the column [1, 0, -1]' 

%%                   
% *Author:* Vladimir Grechka 2012 - 2014
%
% * Test the captions for faults and 3D event-location problems

%%
function [caption, units] = generateCaptions(dim, model, unknowns)
%% Settings 
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

% Unpack 'model' structure
noSou = model.noSou;    anisType = model.anisType;    anisParam = model.anisParam;

if ~(dim == 2  ||  dim == 3)   
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Incorrect dimension of the event-location problem dim = %g \n', dim);
    fprintf('>>> The correct values are dim = 2 or dim = 3 -- STOP \n \n');  
      error(' ');
end;

%% Pieces of captions
indexCij = ['_{11'; '_{12'; '_{13'; '_{14'; '_{15'; '_{16'; ...
                    '_{22'; '_{23'; '_{24'; '_{25'; '_{26'; ...
                            '_{33'; '_{34'; '_{35'; '_{36'; ...
                                    '_{44'; '_{45'; '_{46'; ...
                                            '_{55'; '_{56'; ...
                                                    '_{66'];

%% Set up cell arrays
velISO{1,1} = '{\itV}_{P}';         velISO{2,1} = '{\itV}_{S}';
velAni{1,1} = '{\itV}_{P0}';        velAni{2,1} = '{\itV}_{S0}'; 
aniVTI{1,1} = '\epsilon';           aniVTI{2,1} = '\delta';          aniVTI{3,1} = '\gamma';
aniORT{1,1} = '\epsilon^{(1)}';     aniORT{2,1} = '\epsilon^{(2)}';      
aniORT{3,1} = '\delta^{(1)}';       aniORT{4,1} = '\delta^{(2)}';    aniORT{5,1} = '\delta^{(3)}';
aniORT{6,1} = '\gamma^{(1)}';       aniORT{7,1} = '\gamma^{(2)}';
aniMNC{1,1} = '\zeta^{(1)}';        aniMNC{2,1} = '\zeta^{(2)}';     aniMNC{3,1} = '\zeta^{(3)}';

%% Cij captions
noLayer = size(anisType, 1);
countGlb = 0;
for ilayer=1:noLayer
    ilayerString = strcat('_{ }[', num2str(ilayer), ']}');

    % The numbers of stiffness or anisotropy coefficients by symmetry:
    if strcmp(anisType(ilayer,:), 'ISO') == 1;   noUnknCij = 2;    end;
    if strcmp(anisType(ilayer,:), 'VTI') == 1;   noUnknCij = 5;    end;
    if strcmp(anisType(ilayer,:), 'ORT') == 1;   noUnknCij = 9;    end;
    if strcmp(anisType(ilayer,:), 'MNC') == 1;   noUnknCij = 12;   end;
    if strcmp(anisType(ilayer,:), 'TRI') == 1;   noUnknCij = 21;   end;

    % Derivatives with respect to elastic parameters   
    clear ijUnkn captionCij
    ijUnkn = find(unknowns.CijInd(ilayer,1:noUnknCij) == 0);
    countAdd = length(ijUnkn); 

    if countAdd > 0
        clear captionCij captionAni index
        % Stiffness captions
        if strcmp(anisType(ilayer, :), 'ISO') == 1
            index = indexCij([1,16], :);
        elseif strcmp(anisType(ilayer, :), 'VTI') == 1
            index = indexCij([1,3,12,16,21], :); 
        elseif strcmp(anisType(ilayer, :), 'ORT') == 1
            index = indexCij([1,2,3,7,8,12,16,19,21], :); 
        elseif strcmp(anisType(ilayer, :), 'MNC') == 1
            index = indexCij([1,2,3,6,7,8,11,12,15,16,19,21], :); 
        elseif strcmp(anisType(ilayer, :), 'TRI') == 1
            index = indexCij(1:21, :); 
        else
            fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
            fprintf('>>> Improper value of anisType = %s in layer %g  -- STOP \n \n', ...
                    [anisType(ilayer,:), ilayer]);
              error(' ');
        end;
    
        if strcmp(anisParam(ilayer,:), 'Cij') == 1
            for i=1:size(index,1)
                captionCij{i} = strcat('{\itc}', index(i,:), ilayerString);
            end;
        % Units of dt/dCij
        units(:, (countGlb + 1) : (countGlb + countAdd)) = repmat([3, -2, 0]', 1, countAdd);

        elseif strcmp(anisParam(ilayer,:), 'Ani') == 1
            if strcmp(anisType(ilayer,:), 'ISO') == 1
                captionAni = velISO;
            elseif strcmp(anisType(ilayer,:), 'VTI') == 1
                captionAni = cat(1, velAni, aniVTI);   
            elseif strcmp(anisType(ilayer,:), 'ORT') == 1
                captionAni = cat(1, velAni, aniORT);   
            elseif strcmp(anisType(ilayer,:), 'MNC') == 1
                captionAni = cat(1, velAni, aniORT, aniMNC);   
            else
                fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
                fprintf('>>> Improper value of anisParam = %s \n', anisParam(ilayer,:));
                fprintf('    in layer %g whose symmetry is %s, \n', [ilayer, anisType(ilayer,:)]);
                fprintf('    that is, not ISO or VTI or ORT or MNC \n');
                error('>>> STOP');
            end;

            icount = 0;
            for i=1:size(index,1)
                captionCij{i} = strcat(captionAni(i,:), '_{', ilayerString);
                if i <= length(ijUnkn)
                    if ijUnkn(i) <= 2
                        icount = icount + 1;
                        % Units of dt/d(velocity)
                        units(:, (countGlb + icount)) = [2, -1, 0]';
                    else
                        icount = icount + 1;
                        % Units of dt/d(anisotropy coefficient)
                        units(:, (countGlb + icount)) = [1, 0, 0]';
                    end;         
                end;           
            end;

        else
            fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
            fprintf('>>> Improper value of anisType = %s in layer %g  -- STOP \n \n', ...
                    [anisType(ilayer,:), ilayer,]);
            error(' ');
        end;
                
        [caption{(countGlb + 1) : (countGlb + countAdd)}] = deal(captionCij{ijUnkn});
        countGlb = countGlb + countAdd; 
    end;
    
    % Derivatives with respect to rotation angles   
    if strcmp(anisType(ilayer,:), 'TRI') == 0
        clear ijUnkn captionCij    
        ijUnkn = find(unknowns.CijInd(ilayer, (noUnknCij+1) : (noUnknCij+3)) == 0);
        countAdd = length(ijUnkn); 
        if countAdd > 0
            % Captions for the rotation angles 
            for i=1:3
                captionCij{i} = strcat('\alpha', '^{(', num2str(i), ')}_{', ilayerString);
            end;
            [caption{(countGlb + 1) : (countGlb + countAdd)}] = deal(captionCij{ijUnkn});
            units(:, (countGlb + 1) : (countGlb + countAdd)) = repmat([1, 0, -1]', 1, countAdd);
            countGlb = countGlb + countAdd; 
        end;  
    end;
end;

%% Source-coordinate captions
for isou=1:noSou
    if unknowns.perfInd(isou) == 0
        is = strcat('{', num2str(isou), '}');
        if dim == 2
            caption{countGlb + 1} = strcat('{\itr}_', is);
            caption{countGlb + 2} = strcat('{\itz}_', is);
        else
            caption{countGlb + 1} = strcat('{\xi}_{ 1,', is, '}');
            caption{countGlb + 2} = strcat('{\xi}_{ 2,', is, '}');
            caption{countGlb + 3} = strcat('{\xi}_{ 3,', is, '}');
        end;
        units(:, (countGlb + 1) : (countGlb + dim)) = repmat([1, -1, 0]', 1, dim);
        countGlb = countGlb + dim;
    end;
end;

%% Interface captions
[irow1, icol1] = find(unknowns.interfaceInd == 0);
for iint=1:length(irow1)
    % Determine the interface number
    intNumd = strcat(   '{[', num2str(icol1(iint)), ']}');
    intNumb = strcat('_{ }[', num2str(icol1(iint)), ']}');
    if irow1(iint) == 3
        caption{countGlb + 1} = strcat('{\itd}_', intNumd);
        units(:,countGlb + 1) = [1, -1, 0]';
    elseif irow1(iint) == 4
        caption{countGlb + 1} = strcat('{\itb}_{1', intNumb);
        units(:,countGlb + 1) = [1, 0, -1]';
    elseif irow1(iint) == 5
        caption{countGlb + 1} = strcat('{\itb}_{2', intNumb);  
        units(:,countGlb + 1) = [1, 0, -1]';
    end;
    countGlb = countGlb + 1;
end;

%% Fault captions
[irow2, icol2] = find(unknowns.faultInd == 0);
for iflt=1:length(irow2)
    % Determine the fault number
    intNumd = strcat(   '{[', num2str(icol2(iflt) - 1), ']}');
    intNumb = strcat('_{ }[', num2str(icol2(iflt) - 1), ']}');
    if irow2(iflt) == 3
        caption{countGlb + 1} = strcat('{\itf}_', intNumd);
        units(:,countGlb + 1) = [1, -1, 0]';
    elseif irow2(iflt) == 4
        caption{countGlb + 1} = strcat('{\itn}_{1', intNumb);
        units(:,countGlb + 1) = [1, 0, -1]';
    elseif irow2(iflt) == 5
        caption{countGlb + 1} = strcat('{\itn}_{2', intNumb);  
        units(:,countGlb + 1) = [1, 0, -1]';
    end;
    countGlb = countGlb + 1;
end;

%% Origin-time captions
for isou=1:noSou
    if unknowns.tauInd(isou) == 0
        isouStr = strcat('{', num2str(isou), '}');
        caption{countGlb + 1} = strcat('\tau_', isouStr);
        units(:,countGlb + 1) = [0, 0, 0]';
        countGlb = countGlb + 1;
    end;
end;

end    % of the function