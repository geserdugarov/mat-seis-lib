%% *dFdAniCoef*
% Transform derivatives dF/dCij calculated in TRI media into those with respect to the anisotropy 
% coefficients in VTI, orthorhombic, and monoclinic media  

%% 
% *Input:*

% dCijInp  - [:, 21] matrix containing the rows of derivatives with respect to 21 stiffness
%            coefficients
% Cij      - [6, 6] stiffness matrix in the natural coordinate frame 
% anisType - [1, 3][char] array indicating the layer symmetry:
%            .ISO - isotropy 
%            .VTI - transverse isotropy; the symmetry axis can be tilted
%            .ORT - orthotropy; the symmetry planes can be arbitrarily rotated 
%            .MNC - monoclinic medium with a horizontal symmetry plane; subsequent rotation
%                   is possible
%            .TRI - triclinic symmetry
% indUnkn  - array of the indexes of unknown elastic parameters specified by 'setUnknowns'

%%
% *Output:*

% dAni     - [:, noUnknCij] matrix of the derivatives with respect to the anisotropy 
%            coefficients
% dC       - [21, noUnknCij] matrix of the derivatives dCij/d(Anisotropy coefficients)
%            Parameter 'noUnknCij' is defined in 'reduceDerivMat' for different symmetries

%%
% *Author:* Vladimir Grechka 2012 

%%
% *Comments:* 

% * Equations implemented in this script are derived in 'dFdAniCoef.nb'

% * WARNING: The equations assume strong anisotropy (artFlags.flagWAA = 'SA' in 'art3D') and will 
%            yield incorrect derivatives, possibly causing divergence of the inversion, for weak 
%            anisotropy (artFlags.flagWAA = 'WA' in 'art3D')
%

%%
function [dAni, dC] = dFdAniCoef(dCijInp, Cij, anisType, indUnkn)
%% Settings  
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

% Check the validity of parameter 'anisType' 
if strcmp(anisType, 'ISO') == 0  &&  strcmp(anisType, 'VTI') == 0  &&  ...  
   strcmp(anisType, 'ORT') == 0  &&  strcmp(anisType, 'MNC') == 0  
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Anisotropy parameters are defined for anisType = ISO, VTI, ORT, or MNC \n');
    fprintf('>>> Invalid anisType = %s \n \n', anisType);
      error('>>> STOP');
end;        

%% Compute the most general set of anisotropy coefficients for Cij given in its natural coordinate frame
tmp = cij2grechka(Cij);
Vp0  = tmp(1);   Vs0  = tmp(2); 
eps1 = tmp(3);   eps2 = tmp(4); 
del1 = tmp(5);   del2 = tmp(6);    del3 = tmp(7);
gam1 = tmp(8);   gam2 = tmp(9);
zet1 = tmp(10);  zet2 = tmp(11);   zet3 = tmp(12);

Vp2 = Vp0^2;     Vs2 = Vs0^2;

%% Calculate dC = d Cij / d Ani, see 'dFdAniCoef.nb'        
clear dC sqrt1 sqrt2 expr1 expr2 
if strcmp(anisType, 'ISO') == 1
    % Derivatives for isotropy
    % Parameter vector = [Vp0, Vs0] 
    dC = zeros(21,2);

    % d/dVp0
    dC(1,1)  = 2*Vp0;     dC(2,1)  = dC(1,1);    dC(3,1)  = dC(1,1);
    dC(7,1)  = dC(1,1);   dC(8,1)  = dC(1,1);    dC(12,1) = dC(1,1);

    % d/dVs0
    dC(2,2)  = -4*Vs0;    dC(3,2)  = dC(2,2);    dC(8,2)  = dC(2,2);
    dC(16,2) =  2*Vs0;    dC(19,2) = dC(16,2);   dC(21,2) = dC(16,2);
end;

%%
if strcmp(anisType, 'VTI') == 1
    % Derivatives for VTI
    % Parameter vector = [Vp0, Vs0, eps1, del1, gam1] 
    dC = zeros(21,5);
    sqrt1 = sqrt((Vp0 - Vs0)*(Vp0 + Vs0));
    sqrt2 = sqrt((1 + 2*del1)*Vp2 - Vs2);
    expr1 = 2*(1 + del1)*Vp0*Vs2;

    % d/dVp0
    dC(1,1)  = 2*(1 + 2*eps1)*Vp0;      dC(2,1) = dC(1,1);      dC(7,1) = dC(1,1);
    dC(3,1)  = ((2 + 4*del1)*Vp0^3 - expr1)/(sqrt1*sqrt2);      dC(8,1) = dC(3,1);
    dC(12,1) = 2*Vp0;

    % d/dVs0
    dC(2,2)  = -4*(1 + 2*gam1)*Vs0;
    dC(3,2)  = (-2*Vs0*((1 + del1)*Vp2 - Vs2 + sqrt1*sqrt(Vp2 + 2*del1*Vp2 - Vs2)))/(sqrt1*sqrt2);
    dC(8,2)  = -((Vs0*(sqrt1 + sqrt2)^2)/(sqrt1*sqrt2));
    dC(16,2) = 2*Vs0;
    dC(19,2) = dC(16,2);
    dC(21,2) = 2*(1 + 2*gam1)*Vs0;

    % d/deps1 
    dC(1,3)  = 2*Vp2;   dC(2,3) = dC(1,3);   dC(7,3) = dC(1,3);

    % d/ddel1 
    dC(3,4)  = (Vp2*sqrt1)/sqrt2;   dC(8,4) = dC(3,4);

    % d/dgam1 
    dC(2,5)  = -4*Vs2;
    dC(21,5) = 2*Vs2;
end;

%%
if strcmp(anisType, 'ORT') == 1 || strcmp(anisType, 'MNC') == 1
    % Derivatives for ORT and a portion of derivatives for MNC 
    % Parameter vector: [Vp0, Vs0, eps1, eps2, del1, del2, del3, gam1, gam2] 
    dC = zeros(21,9);
    sqrt1 = sqrt((Vp0 - Vs0)*(Vp0 + Vs0));
    sqrt2 = sqrt((1 + 2*del2)*Vp2 - Vs2);
    expr1 = (1 + 2*gam1)*Vs2;
    expr2 = expr1/(1 + 2*gam2);

    % d/dVp0
    dC(1,1)  = 2*(1 + 2*eps2)*Vp0;
    dC(2,1)  = (2*(1 + 2*eps2)*Vp0*((1 + 2*del3)*(1 + 2*eps2)*Vp2 - (1 + del3)*expr1))/ ...
               (sqrt((1 + 2*eps2)*Vp2 - expr1)*sqrt((1 + 2*del3)*(1 + 2*eps2)*Vp2 - expr1));
    dC(3,1)  = ((2 + 4*del2)*Vp0^3 - 2*(1 + del2)*Vp0*Vs2)/(sqrt1*sqrt2);
    dC(7,1)  = 2*(1 + 2*eps1)*Vp0;
    dC(8,1)  = (2*Vp0*((1 + 2*del1)*(1 + 2*gam2)*Vp2 - (1 + del1)*expr1))/ ...
               ((1 + 2*gam2)*sqrt(Vp2 - expr2)*sqrt((1 + 2*del1)*Vp2 - expr2));
    dC(12,1) = 2*Vp0;

    % d/dVs0
    dC(2,2)  = -(((1 + 2*gam1)*Vs0*(sqrt((1 + 2*eps2)*Vp2 - expr1) + ...
                  sqrt((1 + 2*del3)*(1 + 2*eps2)*Vp2 - expr1))^2)/ ...
                 (sqrt((1 + 2*eps2)*Vp2 - expr1)*sqrt((1 + 2*del3)*(1 + 2*eps2)*Vp2 - expr1)));
    dC(3,2)  = (-2*Vs0*((1 + del2)*Vp2 - Vs2 + sqrt1*sqrt(Vp2 + 2*del2*Vp2 - Vs2)))/(sqrt1*sqrt2);
    dC(8,2)  = -(((1 + 2*gam1)*Vs0*(sqrt(Vp2 - expr2) + sqrt((1 + 2*del1)*Vp2 - expr2))^2)/ ...
                 ((1 + 2*gam2)*sqrt(Vp2 - expr2)*sqrt((1 + 2*del1)*Vp2 - expr2)));
    dC(16,2) = (2*(1 + 2*gam1)*Vs0)/(1 + 2*gam2);
    dC(19,2) = 2*Vs0;
    dC(21,2) = 2*(1 + 2*gam1)*Vs0;

    % d/deps1
    dC(7,3)  = 2*Vp2;

    % d/deps2
    dC(1,4)  = 2*Vp2;
    dC(2,4)  = (2*Vp2*((1 + 2*del3)*(1 + 2*eps2)*Vp2 - (1 + del3)*expr1))/ ...
               (sqrt((1 + 2*eps2)*Vp2 - expr1)*sqrt((1 + 2*del3)*(1 + 2*eps2)*Vp2 - expr1));

    % d/ddel1
    dC(8,5)  = (Vp2*sqrt(Vp2 - expr2))/sqrt((1 + 2*del1)*Vp2 - expr2);

    % d/del2
    dC(3,6)  = (Vp2*sqrt1)/sqrt2;

    % d/del3
    dC(2,7)  = ((1 + 2*eps2)*Vp2*sqrt((1 + 2*eps2)*Vp2 - expr1))/ ...
                    sqrt((1 + 2*del3)*(1 + 2*eps2)*Vp2 - expr1);
    
    % d/dgam1
    dC(2,8)  = -((Vs2*(sqrt((1 + 2*eps2)*Vp2 - expr1) + ...
                  sqrt((1 + 2*del3)*(1 + 2*eps2)*Vp2 - expr1))^2)/ ...
                 (sqrt((1 + 2*eps2)*Vp2 - expr1)*sqrt((1 + 2*del3)*(1 + 2*eps2)*Vp2 - expr1)));
    dC(8,8)  = -((Vs2*(sqrt(Vp2 - expr2) + sqrt((1 + 2*del1)*Vp2 - expr2))^2)/((1 + 2*gam2)* ...
                       sqrt(Vp2 - expr2)*sqrt((1 + 2*del1)*Vp2 - expr2)));
    dC(16,8) = (2*Vs2)/(1 + 2*gam2);
    dC(21,8) = 2*Vs2;

    % d/dgam2
    dC(8,9)  = (expr1*(sqrt(Vp2 - expr2) + sqrt((1 + 2*del1)*Vp2 - expr2))^2)/ ...
               ((1 + 2*gam2)^2*sqrt(Vp2 - expr2)*sqrt((1 + 2*del1)*Vp2 - expr2));
    dC(16,9) = (-2*expr1)/(1 + 2*gam2)^2;
end;

%%
if strcmp(anisType, 'MNC') == 1  
    % Additional MNC derivatives 
    % Parameter vector: [Vp0, Vs0, eps1, eps2, del1, del2, del3, gam1, gam2, zet1, zet2, zet3] 
    dC(:,10:12) = zeros(21,3);

    % d/dVp0
    dC(6,1)  = 2*Vp0*(2*zet1 + zet3);
    dC(11,1) = 2*Vp0*(2*zet2 + zet3);
    dC(15,1) = 2*Vp0*zet3;

    % d/dzet1
    dC(6,10) = 2*Vp2;

    % d/dzet2
    dC(11,11) = 2*Vp2;

    % d/dzet3
    dC(6,12)  = Vp2;   dC(11,12) = dC(6,12);   dC(15,12) = dC(6,12);
end;

%% Derivatives  d / d Ani = (d / d Cij) * (d Cij / d Ani)
dAni = dCijInp*dC(:,indUnkn);
    
end    % of the function
