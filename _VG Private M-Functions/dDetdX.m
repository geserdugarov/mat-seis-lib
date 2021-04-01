%% *dDetdX*
% Differentiate det(F)

%%
% *Input:*

%    F   - square matrix
%   dF   - square matrix containing derivatives of the elements of F

%%
% *Output:*

%   dDet = d det(F) / d X  

%%
% *Author:* Vladimir Grechka 2012 
%
% * Fortran version is published in Obolentseva and Grechka (1989)

%%
function [dDet] = dDetdX(F, dF)
%% Settings  
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

[i1, i2] = size(F);  [j1, j2] = size(dF);

if i1 ~= i2  ||  j1 ~= j2  ||  i1 ~= j1
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('    size(F)  = [%g, %g] \n', [i1, i2]);
    fprintf('    size(dF) = [%g, %g] \n', [j1, j2]);
    fprintf('>>> Either matrix F or dF is not square or their sizes differ \n \n');
      error('>>> STOP');   
end;

%% Compute d[det(F)]/dX 
dDet = 0;
for i=1:size(F,2)
    tmp = F;                  % make matrix tmp equal to F
    tmp(:,i) = dF(:,i);       % replace i-th column of tmp with the corresponding column of dF
    dDet = dDet + det(tmp);
end;

end    % of the function
