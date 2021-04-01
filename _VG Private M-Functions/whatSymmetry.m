%% *whatSymmetry*
% Determine the symmetry of input stiffness matrix Cij 

%%
% *Input:*

% Cij   - [6, 6] stiffness matrix
% tol   - [scalar] tolerance of proximity to a higher-symmetry medium 

%%
% *Output:*

% aniSym - [string] equal to either 'ISO', 'VTI', 'ORT', 'MNC', or 'TRI'
% Rot    - [3, 3] rotation matrix, such that 'CijR = bond(Cij, Rot)' produces 'CijR' in a 
%          crystallographic coordinate frame and 'CijR' is amenable for 'highSymCij' 

%%
% *Author:* Vladimir Grechka 2014

%%
function [aniSym, Rot] = whatSymmetry(Cij, tol)
%% Settings and defaults
% [thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
% narginTrue = nargin(thisFileName);   
% narginchk(narginTrue, narginTrue);

if nargin < 2;   tol = 1.e-4;   end;  % default for the tolerance

%% Determine the medium symmetry
aniSym = 'TRI';

if isISO(Cij, tol) == 1
    aniSym = 'ISO';
    Rot = eye(3);
    return;
else
    
    [flag, Rot] = isTI(Cij, tol);
    if flag == 1
        aniSym = 'VTI';
        return;
    end;

    [flag, Rot] = isORT(Cij, tol);
    if flag == 1
        aniSym = 'ORT';
        return;
    end;

    [flag, Rot] = isMNC(Cij, tol);
    if flag == 1
        aniSym = 'MNC';
        return;
    end;
    
    if strcmp(aniSym, 'TRI') == 1
        Rot = eye(3,3);
        return;
    end;
    
end;

end  % of the function
