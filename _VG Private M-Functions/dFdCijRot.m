%% *dFdCijRot*
% Transform derivatives dF/dCij into those in a rotated coordinate frame 

%%
% *Input:*

% dCijInp - [:, 21] matrix containing the rows of derivatives with respect to 21 stiffness
%           coefficients 
% Rot     - [3, 3] coordinate rotation matrix

%%
% *Output:*

% dCijRot - [:, 21] matrix of the derivatives with respect to elements of the rotated  
%           stiffness matrix given by 'bond(dCijInp, Rot)'

%%
% *Author:* Vladimir Grechka 2012 

%%
function [dCijRot] = dFdCijRot(dCijInp, Rot)
%% Settings  
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

% Check whether rotation matrix Rot is unitary
if isUnitary(Rot)  == 0   
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Non-unitary rotation matrix \n \n');
      error('>>> STOP');   
end;

%% Transformation
dCijRot = NaN(size(dCijInp,1), 21);
for irow=1:size(dCijInp,1)
    dCijVec = dCijInp(irow,:); 

    % Make vector 'dCijVec' a matrix and double its diagonal elements to account for the symmetry
    dCijMat = vec2mat(dCijVec) + diag(diag(vec2mat(dCijVec)));  
               
    % Compute the Bond matrix corresponding to rotation matrix Rot
    [~, B] = bond(dCijMat, Rot, 0); 
    
    % Vector of derivatives with respect to the rotated stiffness components 
    dCijRot(irow,:) = mat2vec(inv(B)*dCijMat*inv(B)').* ...
                      [0.5, 1, 1, 1, 1, 1, 0.5, 1, 1, 1, 1, 0.5, 1, 1, 1, 0.5, 1, 1, 0.5, 1, 0.5];
end;
    
end    % of the function