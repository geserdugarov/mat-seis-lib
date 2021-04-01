%% *dFdRotAngles*
% Find derivatives of the stiffness coefficients CijRot computed as CijRot = bond(Cij, Rot) 
% with respect to the angles rotAngle that define the rotation matrix Rot = loc2glb(rotAngle)
% in the Bond transformation

%%
% *Input:*

% dCijInp    - [:, 21] matrix containing the rows of derivatives with respect to 21 stiffness
%              coefficients
% Cij        - [6, 6] stiffness matrix in the natural coordinate frame 
% rotAngle   - [3, 1] array of three angles specifying the rotation
% indUnknRot - [1, 3] array containing zeros for the elements of 'rotAngle' with respect to 
%              which differentiation is to be performed 

%%
% *Output:*

% dRotAngle - [:, length(rotAngle)] matrix of the derivatives d CijRot/d rotAngle(i)

%%
% *Author:* Vladimir Grechka 2012 

%%
function [dRotAngle] = dFdRotAngles(dCijInp, Cij, rotAngle, indUnknRot)
%% Settings  
[~, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

% Angle increment for numerical differentiation
dangle = 1.e-8;   tol = 1.e-4;   count = 0;

%% Main computation loop
for iangle = indUnknRot
    if iangle == 1  &&  abs(rotAngle(2)) < tol 
        % Avoid division by zero at the singularity of the spherical coordinate system
        count = count + 1;
        dCijRot(:,count) = zeros(21,1);
    else
        if iangle == 1 
            % Account for spherical coordinate system when evaluatng derivative with respect 
            % to azimuth
            da = dangle/abs(sin(rotAngle(2)));
        else
            da = dangle;
        end;
        
        % Approximate derivative with the central difference
        rotAng = rotAngle;
        rotAng(iangle) = rotAng(iangle) + da; 
        if length(rotAng) == 2;   rotAng(3) = 0;   end;
        Rot = loc2glb(rotAng);
        Cup = bond(Cij, Rot);
        
        rotAng = rotAngle;
        rotAng(iangle) = rotAng(iangle) - da; 
        if length(rotAng) == 2;   rotAng(3) = 0;   end;
        Rot = loc2glb(rotAng);
        Cdw = bond(Cij, Rot);
        
        count = count + 1;
        dCijRot(:,count) = (mat2vec(Cup) - mat2vec(Cdw))'/(2*da); 
    end;
end;
         
%% Derivatives with respect to the rotation angles         
dRotAngle = dCijInp*dCijRot;

end    % of the function