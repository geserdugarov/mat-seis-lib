%% *dFdCijSym*
% Transform derivatives dF/dCij calculated in TRI media into those in higher-symmetry media 

%%
% *Input:*

% dCijInp  - [:, 21] matrix containing the rows of derivatives with respect to 21 stiffness
%            coefficients
% anisType - [1, 3][char] array indicating the layer symmetry:
%            .ISO - isotropy 
%            .VTI - transverse isotropy; the symmetry axis can be tilted
%            .ORT - orthotropy; the symmetry planes can be arbitrarily rotated 
%            .MNC - monoclinic medium with a horizontal symmetry plane; subsequent rotation
%                   is possible
%            .TRI - triclinic symmetry
% indUnkn  - array of the indexes of unknown elastic parameters specified in function 
%            'setUnknowns'

%%
% *Output:*

% dCijSym  - [:, length(noUnkn)] matrix of derivatives dF/dCij in high-symmetry media

%%
% *Author:* Vladimir Grechka 2012 - 2014

%%
function [dCijSym] = dFdCijSym(dCijInp, anisType, indUnkn)
%% Settings  
[~, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

clear Mat m2v

%%
if strcmp(anisType, 'ISO') == 1
    % ISO parameters: [c11, c66] and the condition c12 = c11 - 2*c66
    Mat(:,1) = [[ 1  1  1  0  0  0], [ 1  1  0  0  0], [ 1  0  0  0], [ 0  0  0], [ 0  0],  0]';
    Mat(:,2) = [[ 0 -2 -2  0  0  0], [ 0 -2  0  0  0], [ 0  0  0  0], [ 1  0  0], [ 1  0],  1]';
end;
  
%%
if strcmp(anisType, 'VTI') == 1
    % VTI parameters: [c11, c13, c33, c44, c66] and the condition c12 = c11 - 2*c66
    Mat(:,1) = [[ 1  1  0  0  0  0], [ 1  0  0  0  0], [ 0  0  0  0], [ 0  0  0], [ 0  0],  0]';
    Mat(:,2) = [[ 0  0  1  0  0  0], [ 0  1  0  0  0], [ 0  0  0  0], [ 0  0  0], [ 0  0],  0]';
    Mat(:,3) = [[ 0  0  0  0  0  0], [ 0  0  0  0  0], [ 1  0  0  0], [ 0  0  0], [ 0  0],  0]';
    Mat(:,4) = [[ 0  0  0  0  0  0], [ 0  0  0  0  0], [ 0  0  0  0], [ 1  0  0], [ 1  0],  0]';
    Mat(:,5) = [[ 0 -2  0  0  0  0], [ 0  0  0  0  0], [ 0  0  0  0], [ 0  0  0], [ 0  0],  1]';
end;

%%
if strcmp(anisType, 'ORT') == 1
    % ORT parameters: [c11, c12, c13, c22, c23, c33, c44, c55, c66] 
    m2v = [1, 2, 3, 7, 8, 12, 16, 19, 21];
    Mat = full(sparse(m2v, (1:9), ones(1,9), 21, 9));
end;

%%
if strcmp(anisType, 'MNC') == 1
    % MNC parameters: [c11, c12, c13, c16, c22, c23, c26, c33, c36, c44, c55, c66] 
    m2v = [1, 2, 3, 6, 7, 8, 11, 12, 15, 16, 19, 21];
    Mat = full(sparse(m2v, (1:12), ones(1,12), 21, 12));
end;

%%
if strcmp(anisType, 'TRI') == 1
    Mat = diag(ones(1,21), 0);
end;

%% Calculate the Frechet derivatives     
dCijSym = dCijInp*Mat(:,indUnkn);

end    % of the function
