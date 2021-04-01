%% *christoffel*
% Construct tensor G(i,l) = Cij(i,j,k,l)*n1(j)*n2(k), (i, l = 1, 2, 3) that becomes 
% the Christoffel tensor when n1 = n2 is either the wavefront normal or the slowness vector

%%
% *Input:*

% Cij    - [6, 6] stiffness matrix
% n1, n2 - [3, 1] vectors; in general, n1 ~= n2

%%
% *Output:*

% G      - [3, 3] tensor G(i,l) = Cij(i,j,k,l)*n1(j)*n2(k),  (i, l = 1, 2, 3)

%%
% *Author:* Vladimir Grechka 1998 2012
%
% * Fortran version is published in Obolentseva and Grechka (1989)

%%
function G = christoffel(Cij, n1, n2)
%% Settings  
[~, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

flag = isequal(n1, n2);             % check whether n1 = n2
[~, ~, ~, t2m] = indexesVMT(4);     % get the stiffness indexes

%% Construct tensor G 
G = zeros(3, 3);
for i = 1:3
    ind = 1;   
    if flag == 1                % The equality n1 = n2 makes tensor G symmetric and allows 
        ind = i;                % one to reduce the amount of computations
    end;                        
    
    for l = ind:3
        for k = 1:3    
            for j = 1:3   
                G(i,l) = G(i,l) + Cij(t2m(i,j), t2m(k,l))*n1(j)*n2(k);
            end;   
        end;

        if flag == 1  &&  l > i
            G(l,i) = G(i,l);    % symmetrize G
        end;
    end;   
end;

end    % of the function