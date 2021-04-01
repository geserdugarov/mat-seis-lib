%% *scatter2ellipse*
% Approximate a collection of data points with the 95% confidence ellipse or ellipsoid

%%
% *Input:*

% x    - [:, Ndim] array of coordinates of the data points; 
%        Ndim is the dimensionality of the problem  
% flag - [string] equal to either 'eig' (default) or 'svd' that governs the choice of technique
%        used to compute the ellipse
% tol  - [scalar] quantity that stabilizes the subsequent computations by removing near-zero 
%        eigenvalue(s) of the ellipse

%%
% *Output:*

% x0   - [Ndim, 1] coordinates of the center of data   
% A    - [Ndim, Ndim] symmetric positive definite matrix describing the ellipse in the form
%        (x - x0)'*A*(x - x0) = 1
% xvec - [Ndim, Ndim] data eigenvector matrix 
% xval - [Ndim, 1] vector of eigenvalues (or semi-axes) of the ellipse

%%                    
% *Author:* Vladimir Grechka 2012

%%
function [x0, A, xvec, xval] = scatter2ellipse(x, flag, tol)
%% Settings and defaults
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
% narginTrue = nargin(thisFileName);   
% narginchk(narginTrue, narginTrue);

Ndim = size(x,2);

if nargin < 3  ||  isempty(tol)  == 1;    tol = 1.e-6;   end;   
if nargin < 2  ||  isempty(flag) == 1;   flag = 'svd';   end; 

%% Get eigenvalues and eigenvectors of the data 
x0 = (mean(x, 1))';                             % center of gravity of the data
y = x - repmat(x0', size(x,1), 1);              % centered data

if strcmp(flag, 'svd') == 1
    [~, ss, xvec] = svd(y, 0);                  % SVD for computing the eigenvalues and 
    xval = diag(ss)/sqrt(size(x,1));            % principal directions of the data

elseif strcmp(flag, 'eig') == 1
    cv = zeros(Ndim, Ndim);
    for i=1:Ndim;  
        for j=1:i;
            cv(i,j) = dot(y(:,i), y(:,j));      % the covariance matrix for eigenvalue analysis
        if i > j;   cv(j,i) = cv(i,j);   end;
        end;    
    end;
    
    [vec, s1] = eig(cv);                        % the eigenvalues and principal directions 
    [ss, ind] = sort(diag(s1), 'descend');      % of the covariance matrix
    xval = sqrt(ss/size(x,1));
    xvec = vec(:,ind);

else
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Incorrect flag = %s \n', flag);
    fprintf('>>> Legitimate values of variable ''flag'' are ''eig'' or ''svd'' \n \n');
      error('>>> STOP');
end;

%% Construct matrix describing an ellipse
for i=1:size(xval,1)
    if xval(i,1)/max(xval) < tol
       xval(i,1) = tol*max(xval);               % remove near-zero eigenvalues
    end;
end;
A = xvec*(2*diag(xval))^(-2)*xvec';             % matrix that yields (x - x0)'*A*(x - x0) = 1

end   % of the function