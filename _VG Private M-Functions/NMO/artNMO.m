%% *artNMO*
% Computation of the effective NMO cylinder

%% 
% *Input:*

% model.   - input model structure, the same as in 'art3D'
% art.     - output structure from 'art3D'

%%
% *Output:*

% nmo.     - structure containing the following fields:
%   cylInt - [3, 3, 2*noInt+1+noFlt, noTime] array of interval NMO cylinders
%            (*) The third index spans the ray trajectory in reverse order: from a receiver 
%            to a source
%   cylEff - [3, 3, noTime] array of effective NMO cylinders at the source locations

%%
% *Author:* Vladimir Grechka 1998 2014

%%
% *Comments and known issues:*
%
% * Original Matlab version, called 'cyl2ell', is a part of the 'ART' package freely  
%   distributed by the CWP, CSM
%
% * The NMO cylinders for ray trajectories containing fault reflections are set to NaN

%%
function [nmo] = artNMO(model, art)
%% Settings and checks
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

%% Unpack the input structures into local variables
interface = model.interface;    fault  = model.fault;      Cij    = model.Cij.global;
waveCode  = model.waveCode;     xSou   = model.xSou;       xRec   = model.xRec;

%% Array sizes 
noSou   = size(xSou, 2);        noRec   = size(xRec, 2);    noWave  = size(waveCode, 2);
noInt   = size(interface, 2);   noFlt   = size(fault, 2);   noTime  = noSou*noRec*noWave;

if noFlt ~= 0
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Fault reflections are present in the ''art3D'' output \n');
    fprintf('    The NMO cylinders for ray trajectories containing fault reflections \n');
    fprintf('    will be set to NaN -- PAUSE \n \n');  pause;
end;

%% Initialize the output arrays
nmoCylInt0 = NaN(3, 3, 2*noInt+1+noFlt, noWave, noRec, noSou);    % interval NMO cylinders
nmoCylEff0 = NaN(3, 3, noWave, noRec, noSou);                     % effective NMO cylinders
    
%% Computation of the NMO cylinders
%parfor isou = 1:noSou
parfor irec = 1:noRec
    for isou = 1:noSou
        for iwave = 1:noWave
%            [iwave, isou, irec]
            itime = noRec*noWave*(isou - 1) + noWave*(irec - 1) + iwave;  
            
            nmoCylIntLocal = NaN(3, 3, 2*noInt+1+noFlt);
            
            rayCode = art.rcode(:, :, itime); % get the current ray code
            rayCode(:, any(isnan(rayCode), 1)) = [];  % remove columns of NaNs from 'rayCode'
            noSeg = size(rayCode, 2);  % the number of segments for the current ray trajectory
            tSeg  = art.tseg(:, itime);
            time  = 0;

            if isnan(art.time(itime)) == 0  &&  sum(rayCode(4,:)) == 0 
                % Ensure that 
                % (a) the ray trajectory has been computed and 
                % (b) it contains no reflections from faults
                
                %% Loop over ray segments
                for jseg = 1:noSeg
                    iseg = noSeg + 1 - jseg;
                    time = time + tSeg(iseg);  
                    pCurrent = -art.n(:, iseg, itime)/art.Vph(iseg, itime);
                    % Compute the current interval NMO cylinder
                    [Uint, qCurrent] = nmoCylInt(Cij(:, :, rayCode(2, iseg)), ...
                                                 pCurrent, rayCode(1, iseg));
                    nmoCylIntLocal(:, :, jseg) = Uint(:, :);             
                                             
                    if jseg == 1    
                        % The NMO cylinder at the first ray segment
                        Ueff = Uint;

                    else  % for 1 < jseg  &&  jseg <= noSeg
                        % Generalized Dix averaging (Grechka and Tsvankin, 2002, p. 945)
                        % (1) Get the interface normal
                        [~, interfaceNormal] = planeEquation(interface(:, rayCode(3, iseg)));
                        % (2) Compute cross-section of the existing effective NMO cylinder 
                        %     by the next interface
                        Weff = nmoCyl2Ell(Ueff, interfaceNormal);
                        % (3) Get the interval NMO ellipse on the other side of the interface
                        Wint = nmoCyl2Ell(Uint, interfaceNormal);          
                        % (4) Dix-average the NMO ellipses            
                        WeffNew = inv( ...
                            ((time - tSeg(iseg))*inv(Weff) + tSeg(iseg)*inv(Wint)) / time );
                        % (5) Reconstruct the effective NMO cylinder from the new effective ellipse
                        Ueff = nmoEll2Cyl(WeffNew, interfaceNormal, qCurrent);
                    end;
                    
                end;  % of loop over the ray segments
                nmoCylInt0(:, :, :, iwave, irec, isou) = nmoCylIntLocal;             
                nmoCylEff0(:, :, iwave, irec, isou) = Ueff(:, :);

            end;  % of the ' if isnan(art.time(itime)) == 0  &&  sum(rayCode(4,:)) == 0 ' statement  
            
        end;   % of loop over the wave types
    end;
end;

%% Create the output structure
for isou = 1:noSou;
    for irec = 1:noRec;
        for iwave = 1:noWave;
            itime = noRec*noWave*(isou - 1) + noWave*(irec - 1) + iwave; 
            nmo.cylInt(:, :, :, itime) = nmoCylInt0(:, :, :, iwave, irec, isou);
            nmo.cylEff(:, :, itime)    = nmoCylEff0(:, :, iwave, irec, isou);
        end;
    end;
end

end    % of the function
