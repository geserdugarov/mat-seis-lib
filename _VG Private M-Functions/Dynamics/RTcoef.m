%% *RTcoef*
% Calculation of reflection and transmission coefficients (horizontal bounds)

%%
% *Input:*

% CijR        - [6, 6] array of the stiffness matrices of incident wave layer
% CijT        - [6, 6] array of the stiffness matrices of underlying layer
% densR       - [scalar] density of the incident wave layer
% densT       - [scalar] density of the underlying layer
% n           - [3,1] wave normal of the incident wave
% g           - [3,1] normal to the bound
% typeORpolar - [scalar] wave type: 1 - qP, 2 - qS1, 3 - qS2
%               for TI media: 1 - qP, 2 - qSV, 3 - SH
%               OR polarization vector

%%
% *Output:*

% RT     - [1,6] array of the reflection and transmission coefficients
%          [Rp, Rs1, Rs2, Tp, Ts1, Ts2]
% Vph    - [1,6] phase velocities of the reflected and transmitted waves
% U      - [3,6] polarization vectors of the reflected and transmitted waves
% P      - [3,6] slowness vectors of the reflected and transmitted waves

%%
% *Author:* Geser Dugarov 2017

%%
function [RT, Vph, U, P] = RTcoef(CijR, CijT, densR, densT, n, g, typeORpolar)
%% Settings 
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

printErrors = false;

[~, ~, ~, t2m] = indexesVMT(4); % for Cij from [6,6] to [3,3,3,3]
Vph = zeros(1,6);

%% Rotation to y = 0 plane (only horizontal layers)
nProj = [n(1); n(2); 0];
if norm(nProj) > 1.0e-6
    nProj = nProj/norm(nProj);
    if nProj(2) >= 0 % from -180 deg to 180 deg (acos from 0 to 180 only)
        rotAngle = [0; 0; -acos(nProj(1))];
    else
        rotAngle = [0; 0; acos(nProj(1))];
    end
    % Apply rotation  
    if ~isempty(find(rotAngle ~= 0, 1))
        rotMatrix = loc2glb(rotAngle);
        CijR = bond(CijR, rotMatrix);
        CijT = bond(CijT, rotMatrix);
        rotMatrix = loc2glb(-rotAngle);
        n = rotMatrix*n;
        g = rotMatrix*g;
        if numel(typeORpolar) == 3
            typeORpolar = rotMatrix*typeORpolar;
        end
    end
end

%% Processing input data
CijR4 = zeros(3,3,3,3);
CijT4 = zeros(3,3,3,3);
for i = 1:3
    for k = 1:3
        for j = 1:3
            for l = 1:3
                CijR4(i,j,k,l) = CijR(t2m(i,j),t2m(k,l));
                CijT4(i,j,k,l) = CijT(t2m(i,j),t2m(k,l));
            end
        end
    end
end
rhoCijR4 = densR*CijR4;
rhoCijT4 = densT*CijT4;

[V,U] = velPhaseU(CijR,n);
if numel(typeORpolar) == 1
    if isTI(CijR)
        if ~isISO(CijR)
            % need to choose the right wave
            if abs(dot(U(:,2),n)) > 1.0e-6 || abs(dot(U(:,3),n)) > 1.0e-6
                if abs(dot(U(:,2),n)) < abs(dot(U(:,3),n))
                    V = [V(1) V(3) V(2)];
                    U = [U(:,1) U(:,3) U(:,2)];
                end
            else
                temp1 = U(:,2);
                temp2 = U(:,3);
                if sum(abs(temp1) < 1.0e-6) > sum(abs(temp2) < 1.0e-6)
                    V = [V(1) V(3) V(2)];
                    U = [U(:,1) U(:,3) U(:,2)];
                end
            end
        end
        U0 = U(:,typeORpolar); % polarization of the incident wave
        V0 = V(typeORpolar);   % phase velocity of the incident wave
        p0 = n/V0;             % slowness vector of the incident wave
    else
        U0 = U(:,typeORpolar);
        V0 = V(typeORpolar);
        p0 = n/V0;
    end
elseif numel(typeORpolar) == 3 && abs(norm(typeORpolar)-1) < 1.0e-10
    for i = 1:3
        % could be directed in different ways
        if norm(U(:,i) - typeORpolar) < 1.0e-6
            U0 = U(:,i);
            V0 = V(i);
            break;
        elseif  norm(U(:,i) + typeORpolar) < 1.0e-6
            U0 = -U(:,i);
            V0 = V(i);
            break;
        end
    end
    p0 = n/V0;
else
    fprintf('ERROR: Function ''%s''\n',fullfile(thisFolderName, [thisFileName, '.m']));
    error('Wrong format of ''typeORpolar'' input variable.');
end

%% Calculating reflected and transmitted slowness vectors and polarizations
if isISO(CijR)
    Vp = sqrt(CijR(3,3));
    Vs = sqrt(CijR(5,5));
    solR = -sign(n(3)).*sqrt([1/Vp^2-p0(1)^2; 1/Vs^2-p0(1)^2; 1/Vs^2-p0(1)^2]);
else
    if isTI(CijR)
        % need simple version !!!
        syms p3;
        pRT = vpa([p0(1) p0(2) p3]);
        christR = vpa(reshape(pRT*reshape(shiftdim(reshape(pRT*reshape(shiftdim(CijR4,1),3,3*3*3),3,3,3),1),3,3*3),3,3));
        solR = double(vpasolve(det(christR - eye(3)) == 0, p3));
    else
        syms p3;
        % for horizontal planar bound x1 and x2 coordinates of 
        % reflected/transmitted wave slowness similar to the x1 and x2 
        % coordinates of incident wave slowness
        pRT = vpa([p0(1) p0(2) p3]);
        % forming Christoffel matrix with p3 (unknown)
    %     for i = 1:3
    %         for k = 1:3
    %             for j = 1:3
    %                 for l = 1:3
    %                     christR(i,k) = christR(i,k) + CijR4(i,j,k,l)*pRT(j)*pRT(l);
    %                 end
    %             end
    %         end
    %     end
    %     % C_ijkl -> C_jkli
    %     shiftdim(CijT4,1)
    %     % C_jkli -> C_(j)(kli)
    %     reshape(shiftdim(CijT4,1),3,3*3*3)
    %     % X_(kli) = pRT_j*C_(j)(kli)
    %     pRT*reshape(shiftdim(CijT4,1),3,3*3*3)
    %     % X_(kli) -> X_kli
    %     reshape(X_(kli),3,3,3)
    %     % X_kli -> X_lik
    %     shiftdim(reshape(X_(kli),3,3,3),3,3,3),1)
    %     % X_lik -> X_(l)(ik)
    %     reshape(shiftdim(reshape(X_(kli),3,3,3),3,3,3),1),3,3*3)
    %     % Y_(ik) = pRT_l*X_(l)(ik)
    %     pRT*reshape(shiftdim(reshape(X_(kli),3,3,3),3,3,3),1),3,3*3)
    %     % Y_(ik)-> Y_ik
    %     reshape(Y_(ik),3,3)
        christR = vpa(reshape(pRT*reshape(shiftdim(reshape(pRT*reshape(shiftdim(CijR4,1),3,3*3*3),3,3,3),1),3,3*3),3,3));
        % solving equation Cijkl*pj*pl/rho - E == 0, p - slowness vector
        solR = double(vpasolve(det(christR - eye(3)) == 0, p3));
    end
end
if isISO(CijT)
    Vp = sqrt(CijT(3,3));
	Vs = sqrt(CijT(5,5));
    solT = sign(n(3)).*sqrt([1/Vp^2-p0(1)^2; 1/Vs^2-p0(1)^2; 1/Vs^2-p0(1)^2]);
else
    if isTI(CijT)
        % need simple version !!!
        syms p3;
        pRT = vpa([p0(1) p0(2) p3]);
        christT = vpa(reshape(pRT*reshape(shiftdim(reshape(pRT*reshape(shiftdim(CijT4,1),3,3*3*3),3,3,3),1),3,3*3),3,3));
        solT = double(vpasolve(det(christT - eye(3)) == 0, p3));
    else
        syms p3;
        pRT = vpa([p0(1) p0(2) p3]);
        christT = vpa(reshape(pRT*reshape(shiftdim(reshape(pRT*reshape(shiftdim(CijT4,1),3,3*3*3),3,3,3),1),3,3*3),3,3));
        solT = double(vpasolve(det(christT - eye(3)) == 0, p3));
    end
end

% if there are some complex values
if max(abs(imag(solR))) < 1.0e-8
    solR = real(solR);
else
    if printErrors
        fprintf('RTcoef: There are complex values after searching slowness vectors (R).\n');
    end
end
if max(abs(imag(solT))) < 1.0e-8
    solT = real(solT);
else
    if printErrors
        fprintf('RTcoef: There are complex values after searching slowness vectors (T).\n');
    end
end
signSolR = sign(solR);
ind = abs(real(sign(solR))) < 1.0e-8;
signSolR(ind) = imag(signSolR(ind));
signSolT = sign(solT);
ind = abs(real(sign(solT))) < 1.0e-8;
signSolT(ind) = imag(signSolT(ind));

% choosing right directions
if sign(n(3)) > 0
    solR = solR(signSolR < 0);
    solT = solT(signSolT > 0);
else
    solR = solR(signSolR > 0);
    solT = solT(signSolT < 0);
end
if numel(solR) < 3
    fprintf('ERROR: Function ''%s''\n',fullfile(thisFolderName, [thisFileName, '.m']));
    error('Can not find Christoffel matrix eigenvalues for reflected waves.\n');
end
if numel(solT) < 3
    fprintf('ERROR: Function ''%s''\n',fullfile(thisFolderName, [thisFileName, '.m']));
    error('Can not find Christoffel matrix eigenvalues for transmitted waves.\n');
end

% sorting values: [P-wave, S1-wave, S2-wave]
[~,I] = max([sum(abs(solR - solR(1))) sum(abs(solR - solR(2))) sum(abs(solR - solR(3)))]); % index for P-wave
if abs(solR(mod(I,3)+1)) < abs(solR(mod(I+1,3)+1)) % other two indexes
    solR = [solR(I) solR(mod(I,3)+1)   solR(mod(I+1,3)+1)];
else
    solR = [solR(I) solR(mod(I+1,3)+1) solR(mod(I,3)+1)];
end
[~,I] = max([sum(abs(solT - solT(1))) sum(abs(solT - solT(2))) sum(abs(solT - solT(3)))]); % index for P-wave
if abs(solT(mod(I,3)+1)) < abs(solT(mod(I+1,3)+1)) % other two indexes
    solT = [solT(I) solT(mod(I,3)+1)   solT(mod(I+1,3)+1)];
else
    solT = [solT(I) solT(mod(I+1,3)+1) solT(mod(I,3)+1)];
end

% forming slowness vectors
pRp  = [p0(1); p0(2); solR(1)];
pRs1 = [p0(1); p0(2); solR(2)];
pRs2 = [p0(1); p0(2); solR(3)];
pTp  = [p0(1); p0(2); solT(1)];
pTs1 = [p0(1); p0(2); solT(2)];
pTs2 = [p0(1); p0(2); solT(3)];
if abs(solR(2)-solR(3)) < 1.0e-15
    pRs2 = pRs1;
end
if abs(solT(2)-solT(3)) < 1.0e-15
    pTs2 = pTs1;
end
nRp  = pRp/norm(pRp);
nRs1 = pRs1/norm(pRs1);
nRs2 = pRs2/norm(pRs2);
nTp  = pTp/norm(pTp);
nTs1 = pTs1/norm(pTs1);
nTs2 = pTs2/norm(pTs2);

% forming polarization vectors
if isISO(CijR)
    Vp = sqrt(CijR(3,3));
    Vs = sqrt(CijR(5,5));
    UrP  = polarVect(CijR,pRp);
    UrP2 = polarVect(CijR,pRs1*Vs/Vp); %for S-wave there is another P-wave polarization
    if norm(nRs2-g) < 1.0e-6 || norm(nRs2+g) < 1.0e-6
        [~, U] = velPhaseU(CijR,nRs1);
        UrS1 = U(:,2);
        UrS2 = U(:,3);
    else
        UrS2 = cross(pRs2,g);
        UrS2 = UrS2/norm(UrS2);
        UrS1 = cross(UrP2,UrS2);
        UrS1 = UrS1/norm(UrS1);
    end
    % not always right hand system
    if acos(dot(cross( UrS1,UrS2 ),UrP)) > pi/2
        UrP  = -UrP;
    end
    clear UrP2;
else
    UrP  = polarVect(CijR,pRp);
    UrS1 = polarVect(CijR,pRs1);
    if size(UrS1,2) > 1
        % singularity point, any 2 perpendicular vectors are suit
        [~, U] = velPhaseU(CijR,nRs1);
        if norm(U(:,1)-g) < 1.0e-6 || norm(U(:,1)+g) < 1.0e-6
            UrS1 = U(:,2);
            UrS2 = U(:,3);
        else
            UrS2 = cross(U(:,1),g);
            UrS2 = UrS2/norm(UrS2);
            UrS1 = cross(U(:,1),UrS2);
            UrS1 = UrS1/norm(UrS1);
        end
    else
        UrS2 = polarVect(CijR,pRs2);
    end
end
if isISO(CijT)
    Vp = sqrt(CijT(3,3));
    Vs = sqrt(CijT(5,5));
    UtP  = polarVect(CijT,pTp);
    UtP2 = polarVect(CijT,pTs1*Vs/Vp); %for S-wave there is another P-wave polarization
    if norm(nTs2-g) < 1.0e-6 || norm(nTs2+g) < 1.0e-6
        [~, U] = velPhaseU(CijT,nTs1);
        UtS1 = U(:,2);
        UtS2 = U(:,3);
    else
        UtS2 = cross(pTs2,g);
        UtS2 = UtS2/norm(UtS2);
        UtS1 = cross(UtP2,UtS2);
        UtS1 = UtS1/norm(UtS1);
    end
    % not always right hand system
    if acos(dot(cross( UtS1,UtS2 ),UtP)) > pi/2
        UtP  = -UtP;
    end
    clear UtP2;
else
    UtP  = polarVect(CijT,pTp);
    UtS1 = polarVect(CijT,pTs1);
    if size(UtS1,2) > 1
        % singularity point, any 2 perpendicular vectors are suit
        [~, U] = velPhaseU(CijT,nTs1);
        if norm(U(:,1)-g) < 1.0e-6 || norm(U(:,1)+g) < 1.0e-6
            UtS1 = U(:,2);
            UtS2 = U(:,3);
        else
            UtS2 = cross(U(:,1),g);
            UtS2 = UtS2/norm(UtS2);
            UtS1 = cross(U(:,1),UtS2);
            UtS1 = UtS1/norm(UtS1);
        end
    else
        UtS2 = polarVect(CijT,pTs2);
    end
end

% checking norm of polarization vectors
% for complex case (U,U)=1, not (U,U*)=1 
UrP  = UrP  / sqrt(UrP(1)^2  + UrP(2)^2  + UrP(3)^2);
UrS1 = UrS1 / sqrt(UrS1(1)^2 + UrS1(2)^2 + UrS1(3)^2);
UrS2 = UrS2 / sqrt(UrS2(1)^2 + UrS2(2)^2 + UrS2(3)^2);
UtP  = UtP  / sqrt(UtP(1)^2  + UtP(2)^2  + UtP(3)^2);
UtS1 = UtS1 / sqrt(UtS1(1)^2 + UtS1(2)^2 + UtS1(3)^2);
UtS2 = UtS2 / sqrt(UtS2(1)^2 + UtS2(2)^2 + UtS2(3)^2);

% checking polarization directions
% two vector directions changing for saving right-handed coordinate system
if acos(dot( UrP,nRp )) > pi/2
    UrP  = -UrP;
    UrS1 = -UrS1;
end
if acos(dot( UtP,nTp )) > pi/2
    UtP  = -UtP;
    UtS1 = -UtS1;
end
if abs(acos(dot( UrS1,nRs1 ))) > abs(acos(dot( UrS2,nRs2 )))
    if acos(dot( UrS1,nRs1 )) > pi/2
        UrS1 = -UrS1;
        UrS2 = -UrS2;
    end
else
    if acos(dot( UrS2,nRs2 )) > pi/2
        UrS1 = -UrS1;
        UrS2 = -UrS2;
    end
end
if abs(acos(dot( UtS1,nTs1 ))) > abs(acos(dot( UtS2,nTs2 )))
    if acos(dot( UtS1,nTs1 )) > pi/2
        UtS1 = -UtS1;
        UtS2 = -UtS2;
    end
else
    if acos(dot( UtS2,nTs2 )) > pi/2
        UtS1 = -UtS1;
        UtS2 = -UtS2;
    end
end
if acos(dot( UrS1,U0 )) > pi/2
    UrS1 = -UrS1;
    UrS2 = -UrS2;
end
if acos(dot( UtS1,U0 )) > pi/2
    UtS1 = -UtS1;
    UtS2 = -UtS2;
end

%% Calculating reflection and transmission coefficients
% CpuTp  = zeros(3,1); CpuRp  = zeros(3,1); CpuTs1 = zeros(3,1); CpuRs1 = zeros(3,1);
% CpuTs2 = zeros(3,1); CpuRs2 = zeros(3,1); Cpu0   = zeros(3,1);
% for i = 1:3
%     for j = 1:3
%         for k = 1:3
%             CpuTp(i)  = CpuTp(i)  + rhoCijT(t2m(i,3), t2m(j,k))*pTp(j)*UtP(k);
%             CpuRp(i)  = CpuRp(i)  + rhoCijR(t2m(i,3), t2m(j,k))*pRp(j)*UrP(k);
%             CpuTs1(i) = CpuTs1(i) + rhoCijT(t2m(i,3), t2m(j,k))*pTs1(j)*UtS1(k);
%             CpuRs1(i) = CpuRs1(i) + rhoCijR(t2m(i,3), t2m(j,k))*pRs1(j)*UrS1(k);
%             CpuTs2(i) = CpuTs2(i) + rhoCijT(t2m(i,3), t2m(j,k))*pTs2(j)*UtS2(k);
%             CpuRs2(i) = CpuRs2(i) + rhoCijR(t2m(i,3), t2m(j,k))*pRs2(j)*UrS2(k);
%             Cpu0(i)   = Cpu0(i)   + rhoCijR(t2m(i,3), t2m(j,k))*p0(j)*U0(k);
%         end
%     end
% end
CpuTp  =  UtP.'*reshape( pTp.'*reshape([0 0 1]*reshape(shiftdim(rhoCijT4,1),3,3*3*3),3,3*3),3,3);
CpuRp  =  UrP.'*reshape( pRp.'*reshape([0 0 1]*reshape(shiftdim(rhoCijR4,1),3,3*3*3),3,3*3),3,3);
CpuTs1 = UtS1.'*reshape(pTs1.'*reshape([0 0 1]*reshape(shiftdim(rhoCijT4,1),3,3*3*3),3,3*3),3,3);
CpuRs1 = UrS1.'*reshape(pRs1.'*reshape([0 0 1]*reshape(shiftdim(rhoCijR4,1),3,3*3*3),3,3*3),3,3);
CpuTs2 = UtS2.'*reshape(pTs2.'*reshape([0 0 1]*reshape(shiftdim(rhoCijT4,1),3,3*3*3),3,3*3),3,3);
CpuRs2 = UrS2.'*reshape(pRs2.'*reshape([0 0 1]*reshape(shiftdim(rhoCijR4,1),3,3*3*3),3,3*3),3,3);
Cpu0   =   U0.'*reshape(  p0.'*reshape([0 0 1]*reshape(shiftdim(rhoCijR4,1),3,3*3*3),3,3*3),3,3);

% equation system
% UtP(1)*Tp + UtS1(1)*Ts1 + UtS2(1)*Ts2 == U0(1) + UrP(1)*Rp + UrS1(1)*Rs1 + UrS2(1)*Rs2;
% UtP(2)*Tp + UtS1(2)*Ts1 + UtS2(2)*Ts2 == U0(2) + UrP(2)*Rp + UrS1(2)*Rs1 + UrS2(2)*Rs2;
% UtP(3)*Tp + UtS1(3)*Ts1 + UtS2(3)*Ts2 == U0(3) + UrP(3)*Rp + UrS1(3)*Rs1 + UrS2(3)*Rs2;
% CpuTp(1)*Tp + CpuTs1(1)*Ts1 + CpuTs2(1)*Ts2 == Cpu0(1) + CpuRp(1)*Rp + CpuRs1(1)*Rs1 + CpuRs2(1)*Rs2;
% CpuTp(2)*Tp + CpuTs1(2)*Ts1 + CpuTs2(2)*Ts2 == Cpu0(2) + CpuRp(2)*Rp + CpuRs1(2)*Rs1 + CpuRs2(2)*Rs2;
% CpuTp(3)*Tp + CpuTs1(3)*Ts1 + CpuTs2(3)*Ts2 == Cpu0(3) + CpuRp(3)*Rp + CpuRs1(3)*Rs1 + CpuRs2(3)*Rs2;

% solving of the system
A = [-UrP(1)   -UrS1(1)   -UrS2(1)   UtP(1)   UtS1(1)   UtS2(1)  ;
     -UrP(2)   -UrS1(2)   -UrS2(2)   UtP(2)   UtS1(2)   UtS2(2)  ;
     -UrP(3)   -UrS1(3)   -UrS2(3)   UtP(3)   UtS1(3)   UtS2(3)  ;
     -CpuRp(1) -CpuRs1(1) -CpuRs2(1) CpuTp(1) CpuTs1(1) CpuTs2(1);
     -CpuRp(2) -CpuRs1(2) -CpuRs2(2) CpuTp(2) CpuTs1(2) CpuTs2(2);
     -CpuRp(3) -CpuRs1(3) -CpuRs2(3) CpuTp(3) CpuTs1(3) CpuTs2(3)];
b = [U0(1); U0(2); U0(3); Cpu0(1); Cpu0(2); Cpu0(3)];
RT = (A\b).';

%% Processing output data
U  = [UrP UrS1 UrS2 UtP UtS1 UtS2];
P  = [pRp pRs1 pRs2 pTp pTs1 pTs2];

end    % of the function

%%