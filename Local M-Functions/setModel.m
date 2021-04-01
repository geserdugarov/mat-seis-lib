%% *setModel*
% Set up a model for ray-tracing calculations

%%
% *Input:*

% anisType            - [:, 3][char] character array specifying the layer symmetry:
%                       . ISO - isotropy, 
%                       . VTI - transverse isotropy (the symmetry axis is not necessarily vertical), 
%                       . ORT - orthotropy,
%                       . MNC - monoclinic symmetry, and
%                       . TRI - triclinic symmetry

%% 
% *Output:*

% model.              -  structure containing the following fields:
%
%     Cij.global      - [6, 6, noLayer] array of the stiffness matrices in global coordinates
%     Cij.symmetry    - [6, 6, noLayer] array of the stiffness matrices in the symmetry (principal) axes
%     rotation.matrix - [3, noLayer] array of three angles specifying the rotation of of the
%                       stiffness tensors from 'Cij.symmetry' to 'Cij.global'
%     rotation.matrix - [3, 3, noLayer] array of rotation matrices given by
%                       rotation.matrix(:,:,ilayer) = loc2glb(rotation.matrix(:,ilayer))  
%     density         - [1, noLayer] array of the densities
%     noLayer         = size(anisType, 1) - the number of model layers 
%
%     interface       - [5, noInt] interface array  
%     noInt           - [scalar] number of interfaces
%
%     fault           - [5, noFlt] fault array  
%     noFlt           - [scalar] number of faults
%
%     xSou            - [3, noSou] array of the source coordinates
%     noSou           - [scalar] number of sources
%     xRec            - [3, noRec] array of the receiver coordinates
%     noRec           - [scalar] number of receivers
%     xWell           - [2, 1] array of lateral coordinates of a vertical observation well
%     noPerf          - [scalar] number of perforation shots
%     noEvnt          - [scalar] number of microseismic events
%                       noSou = noPerf + noEvnt
%     indexPerf       - [1, noSou] array of the indexes of perforation shots, which are found as 
%                       find(indexPerf == 1); other elements of indexPerf are 0
%     indexEvnt       - [1, noSou] array of the indexes of microseismic events, which are found as 
%                       find(indexEvnt == 1); other elements of indexEvnt are 0
%
%     tau             - [1, noSou] array of the origin times of perforation shots and microseismic events
%
%     waveCode        - [4, noWave] wave-code array
%     noWave          - [scalar] number of wave codes

%%
% *Author:* Vladimir Grechka 2012

%%
function [model] = setModel(anisType)
%% Settings 
[~, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

model.anisType = anisType;

%% Input Cij's 
[model.Cij.global, model.Cij.symmetry, model.rotation.angles, model.rotation.matrix, ...
    model.density, model.noLayer] = setCij(anisType);    

%% Input model interfaces and faults
[model.interface, model.noInt] = setInterface;   
[model.fault,     model.noFlt] = setFault;   

%% Input positions of sources and receivers
[model.xSou, model.noSou, model.xRec, model.noRec, model.xWell, ...
    model.noPerf, model.noEvnt, model.indexPerf, model.indexEvnt] = setSouRec;   

%% Input origin the times
[model.tau] = setTau(model.noSou);    

%% Input the wave codes
[model.waveCode, model.noWave] = setWaveCode;    

end    % of the function
