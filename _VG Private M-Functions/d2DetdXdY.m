%% *d2DetdXdY*
% Compute second-order partial derivative of det(F)

%%
% *Input:*

%   F    - square matrix
%   FX   - square matrix of derivatives dF / dX
%   FX   - square matrix of derivatives dF / dY
%   FXY  - square matrix of derivatives d^2 F / (dX dY)

%%
% *Output:*

%   d2Det = d^2 det(F) / (dX dY) 

%%
% *Author:* Vladimir Grechka 2014 
%
% * Fortran version is published in Obolentseva and Grechka (1989)

%%
function [d2Det] = d2DetdXdY(F, FX, FY, FXY) 
%% Settings  
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

[i1, i2] = size(F);  [ix1, ix2] = size(FX);  [iy1, iy2] = size(FY);  [ixy1, ixy2] = size(FXY);
ind = [i1, i2, ix1, ix2, iy1, iy2, ixy1, ixy2];

if isempty(find(diff(ind) ~= 0, 1)) == 0
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('    size(F)   = [%g, %g] \n', [i1, i2]);
    fprintf('    size(FX)  = [%g, %g] \n', [ix1, ix2]);
    fprintf('    size(FY)  = [%g, %g] \n', [iy1, iy2]);
    fprintf('    size(FXY) = [%g, %g] \n', [ixy1, ixy2]);
    fprintf('>>> Matrices above are either not square or have different sizes \n \n');
      error('>>> STOP');    
end;

%% Compute d^2 det(F) / (dX dY) 
d2Det = 0;
for i=1:size(F,2)
    for j=1:size(F,2)
        M = F;                  % make matrix M equal to F
        if i == j
            M(:,i) = FXY(:,i);  % replace i = j-th column of M with the corresponding column of FXY
        else
            M(:,i) = FX(:,i);   % replace i-th column of M with the corresponding column of FX
            M(:,j) = FY(:,j);   % replace j-th column of M with the corresponding column of FY
        end;
        d2Det = d2Det + det(M);
    end;
end;
        
end    % of the function
