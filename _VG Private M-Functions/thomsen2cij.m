%% *thomsen2cij* 
% Calculate stiffness matrix in a VTI medium from its Thomsen's parameters 

%%
% *Input:*

% Ani - [1, 5] vector of Thomsen anisotropy parameters arranged as
%       [Vp0, Vs0, epsilon, delta, gamma]
% msg - [scalar] turning on (msg = 1, default) and off (msg = 0) a verbose warning message

%%
% *Output:*

% Cij - [6, 6] VTI stiffness matrix

%%
% *Author:*  Vladimir Grechka 1996 2012 - 2014

%%
function [Cij] = thomsen2cij(Ani, msg)
%% Settings and defaults
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
%narginTrue = nargin(thisFileName);   
%narginchk(narginTrue, narginTrue); 

if nargin < 2  ||  msg ~= 0
    msg = 1;   % default for printing verbose messages
end

%% Calculate the stiffnesses
c33 = Ani(1)^2;   c11 = c33*(1 + 2*Ani(3));   
c55 = Ani(2)^2;   c66 = c55*(1 + 2*Ani(5));

cc = (c33 - c55)*(2*Ani(4)*c33 + c33 - c55);
if cc > 0  
    if c33 > c55
        c13 = sqrt(cc) - c55;
    else
        c13 = -sqrt(cc) - c55;
        if msg == 1
            fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
            fprintf('>>> Untypical parameters encountered: c33 (= %g) < c55 (= %g) \n', c33, c55);
        end
    end
else 
    c13 = 0;
    if msg == 1
        fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
        fprintf('>>> Coefficient c13 is set to zero to avoid a complex-valued stiffness matrix \n');
    end
end

%% Construct Cij's 
Cij = zeros(6, 6);
Cij(1,1) = c11;   Cij(1,2) = c11 - 2*c66;   Cij(1,3) = c13;
                  Cij(2,2) = c11;           Cij(2,3) = c13;
                                            Cij(3,3) = c33;
Cij(4,4) = c55;   Cij(5,5) = c55;           Cij(6,6) = c66;
Cij = symMat(Cij);   % fill the symmetric part

end   % of the function