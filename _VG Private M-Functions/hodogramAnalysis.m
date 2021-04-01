%% *hodogramAnalysis*
% Hodogram analysis of 3C seismic data 

%%
% *Input:*

% traceData        - [:, 3] array containing the input 3C trace
% angleUnits       - [char] output angle units, 'rad' or 'deg', for the hodogram fields 
%                    'polAngle', 'polStd', 'azmAngle', and 'azmStd' 

%%
% *Output:*

% hodogram.        - structure containing the following fields:        
%     eigVal       - [3, 1] array of eigenvalues of the particle motion 
%     eigVec       - [3, 3] array of eigenvectors of the particle motion 
%     polAngle     - [scalar] polar angle (rad) of the major eigenvector of 'traceData'
%     polStd       - [scalar] standard deviation (rad) of 'polAngle'
%     azmAngle     - [scalar] azimuth (rad) of the major eigenvector of 'traceData'
%     azmStd       - [scalar] standard deviation (rad) of 'azmAngle'
%     linearity    - [scalar] linearity measure of the particle motion 
%     planarity    - [scalar] planarity measure of the particle motion 
%     amplitudeMax - [scalar] maximum amplitude of the particle motion 
%     amplitudeRMS - [scalar] RMS amplitude of the particle motion 

%%
% *Author:* Vladimir Grechka 2012 - 2014

%%
function [hodogram] = hodogramAnalysis(traceData, angleUnits)
%% Settings and defaults
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

if isempty(traceData) == 1
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Empty trace array \n');    
    fprintf('>>> PAUSE -- Continue? \n');    
    return;
end;

if size(traceData, 2) ~= 3
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> size(traceData, 2) should be equal to 3 rather than %g \n \n', ...
            size(traceData, 2));    
      error('>>> STOP');    
end;

%% Eigenvalue-eigenvector analysis of the data
[x0, A, eigVec, eigVal] = scatter2ellipse(traceData, [], []);

%if max(abs(eigVec(:,1))) ~= max(eigVec(:,1))            % make the maximum component of the major 
%    eigVec(:,1) = -eigVec(:,1);                         % eigenvector positive
%end;

if eigVec(3,1) == abs(eigVec(3,1))                       % make the z-component of the major 
    eigVec(:,1) = -eigVec(:,1);                          % eigenvector negative
end;

if dot(cross(eigVec(:,1), eigVec(:,2)), eigVec(:,3)) < 0
    eigVec(:,3) = -eigVec(:,3);                         % make the eigenvector matrix left-handed
end;

%% Amplitudes
traceEig = dot(traceData, repmat(eigVec(:,1)', size(traceData, 1), 1), 2);
amplitudeMax = max(abs(traceEig));                                  % maximum amplitude
amplitudeRMS = norm(traceEig)/sqrt(size(traceData, 1));             % RMS amplitude 

%% The linearity and planarity measures
% Grechka and Mateeva (2007, Geophysics)
linearity = (3*max(eigVal)/sum(eigVal) - 1)/2;         

% Moriya and Niitsuma (1996, Geophysics)
% linearity = ((eigVal(1) - eigVal(2))^2 + ...           
%              (eigVal(1) - eigVal(3))^2 + ...  
%              (eigVal(2) - eigVal(3))^2)/(2*sum(eigVal)^2);

% Jurkevics (1988, Polarization analysis of 3C array data: BSSA, 78, no.5, 1725-1743)
% linearity = 1 - (eigVal(2) + eigVal(3))/(2*eigVal(1));  
planarity = 1 - 2*eigVal(3)/(eigVal(1) + eigVal(2));       

% Note: The linearities defined by Grechka and Mateeva (2007) and Moriya and Niitsuma (1996) are  
%       guaranteed to be no greater than the planarity (see 'A_LinPlanTest.m'), whereas the 
%       linearity defined by Jurkevics (1988) can be (unphysically) greater than the planarity

%% Polar and azimuthal angles
polAngle = atan2(sqrt(eigVec(1,1)^2 + eigVec(2,1)^2), eigVec(3,1));  % (De Meersman et al., 2006,
azmAngle = atan2(eigVec(2,1), eigVec(1,1));                          % BSSA, 96, No 6, 2415-2430)

%% Standard deviations in the polar and azimuthal angles
rot = [cos(azmAngle), -sin(azmAngle), 0; ...
       sin(azmAngle),  cos(azmAngle), 0; 0, 0, 1];
traceDataXZ = traceData*rot;                            % rotate data to the [x, z]-plane

[uDataXZ, sDataXZ, vDataXZ] = svd([traceDataXZ(:,1), traceDataXZ(:,3)], 0);
polStd = atan2(sDataXZ(2,2), sDataXZ(1,1));

[uDataXY, sDataXY, vDataXY] = svd([traceDataXZ(:,1), traceDataXZ(:,2)], 0); 
azmStd = atan2(sDataXY(2,2), sDataXY(1,1));

%% Assign the fields of output structure
hodogram.center   = x0;      
hodogram.matrix   = A;      
hodogram.eigVal   = eigVal;      
hodogram.eigVec   = eigVec;    
hodogram.polAngle = u2u(polAngle, 'rad', angleUnits);    
hodogram.polStd   = u2u(polStd,   'rad', angleUnits);
hodogram.azmAngle = u2u(azmAngle, 'rad', angleUnits);    
hodogram.azmStd   = u2u(azmStd,   'rad', angleUnits);
hodogram.linearity    = linearity;   
hodogram.planarity    = planarity;   
hodogram.amplitudeMax = amplitudeMax;    
hodogram.amplitudeRMS = amplitudeRMS;

end    % of the function