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
% *Author:* Geser Dugarov 2018

%%
function G = christoffel(Cij, n1, n2)

[~, ~, ~, t2m] = indexesVMT(4);     % get the stiffness indexes

if size(n1,1) == 1
    n1 = n1';
end
if size(n2,1) == 1
    n2 = n2';
end

Cij4 = zeros(3,3,3,3);
for i = 1:3
    for j = 1:3
        for k = 1:3
            for l = 1:3
                Cij4(i,j,k,l) = Cij(t2m(i,j),t2m(k,l));
            end
        end
    end
end

%% Construct tensor G 
G = double(vpa(reshape(n1'*reshape(shiftdim(reshape(n1'*reshape(shiftdim(Cij4,1),3,3*3*3),3,3,3),1),3,3*3),3,3)));

end    % of the function