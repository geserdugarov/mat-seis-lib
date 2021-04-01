%% *frechetSVD*
% SVD analysis of Frechet-derivative matrix to determine the scaling of the units of time, length,
% and angle that minimizes its condition number

%%
% *Input:*

% fig           - [scalar] number of an existing figure
% Mat           - [:, :] the Frechet-derivative matrix whose SVD analysis is to be performed
% units         - [3, :] matrix obtained from function 'generateCaptions'
% svdFlag       - [4, 1] array containing ones and zeros to indicate which of the four
%                 outputs of the SVD analysis (the ones) is to be saved
% figPosition   - [4, 1] vector specified the figure position and size
% figRHS        - [cell] array of captions constructed in 'generateCaptions'
% figNamePrefix - [char] prefix of the file name for saving generated figures
% saveFigFlag   - [3, 1] array containing ones and zeros used in 'saveFigure'
% plotEigVec    - [scalar] number determining the maximum size of eigenvector matrix 
%                 to be plotted
% sizeSV        - [scalar] size of circles indicating the singular values
% widthSV       - [scalar] width of line connecting the singular values
% verboseFlag   - [scalar] flag, equal to 0 or 1, to indicate verbose print during optimization

%%
% *Output:*

% fignext       - [scalar] number of the next figure
% scaleTLA      - [3, 1] array of computed time, length, and angle scalings in the units of
%                 quantities in 'Mat'
% MatOut        - [size(Mat)] scaled version of input matrix 'Mat'
% scaleOut      - [1, size(Mat,2)] array of the inverse scaling factors for the parameter vector
% numCondOut    - [scalar] condition number of the properly scaled matrix 'Mat'

%%
% *Author:* Vladimir Grechka 2012 - 2014

%%
% *2012.12.19 Note:* The column normalization suggested in the PCH book is much faster than 
% the optimization for the TLA scaling factors and often numerically better, that is, it yields a
% smalle condition number. Perhaps I should stick to it and remove 'condNumOptim.m' once it again 
% produces NaNs or Infs 

% *2013.09.11 Warning Note:* If a column of input matrix is numerically zero up to round-off 
% errors, the PCH scaling amplifies it with a very large scaling coefficient and may lead to 
% erroneous conclusions of the SVD analysis. Temporarily, I deal with the issue by introducing     
% the parameter 'tolMax'.

%%
function [fignext, scaleTLA, MatOut, scaleOut, numCondOut] = ...
      frechetSVD(fig, Mat, units, svdFlag, figPosition, figRHS, figNamePrefix, ...
                 saveFigFlag, plotEigVec, sizeSV, widthSV, verboseFlag)

global scaleMat scaleCond numCond

%% Settings and defaults
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

scaleTLA = [1 1 1];    MatOut = Mat;    scaleOut = ones(1,size(Mat,2));

if isempty(plotEigVec)  == 1;   plotEigVec = 50;   end;
if isempty(verboseFlag) == 1;   verboseFlag = 0;   end;

figXlabel = 'Parameter_{ } combination';
figYlabel = 'Log_{ 10} of normalized singular value';

%% SVD before scaling
if svdFlag(1) == 1  ||  svdFlag(2) == 1
    [~, s, w] = svd(Mat,0);  
    
    if svdFlag(1) == 1
        permuteOrder = 1:size(Mat,2);
        fig = fig + 1;
        setFigure(fig, 2, figPosition, {figNamePrefix; 'SVD before scaling'}, ...
                  figXlabel, figYlabel, [], [], []);
        plotSVD(fig, diag(s), w, figRHS, permuteOrder, plotEigVec, [], [], [], sizeSV, widthSV);
        saveFigure(fig, fullfile('Figures', [figNamePrefix, 'svd-01']), saveFigFlag);
    end;
    
    if svdFlag(2) == 1
        permuteOrder = enhanceDiagDom(w);
        [~, s, w] = svd(Mat(:,permuteOrder),0);
        fig = fig + 1;
        setFigure(fig, 2, figPosition, {figNamePrefix; 'Sorted SVD before scaling'}, ...
                  figXlabel, figYlabel, [], [], []);
        plotSVD(fig, diag(s), w, figRHS, permuteOrder, plotEigVec, [], [], [], sizeSV, widthSV);
        saveFigure(fig, fullfile('Figures', [figNamePrefix, 'svd-02']), saveFigFlag);
    end;
end;

%% SVD after scaling
if svdFlag(3) == 1  ||  svdFlag(4) == 1
    % Adjust the units to get the smallest condition number
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('\n>>> Minimizing the condition number... \n');
    scaleTLA = fmincon('condNumOptim', [1 1 1], [], [], [], [], ...
                        1.e-6*[1 1 1], 1.e+6*[1 1 1], [], [], Mat, units, verboseFlag);
    numCondOut = numCond;

    % Scaling proposed in Hansen et al. (2012) 
    [MatOut, scaleOut, scaleFlag, numCondPCH] = condNumHansen(Mat, 0);
        
    if numCondPCH < numCondOut  &&  scaleFlag == 1
        numCondOut = numCondPCH;
        scaleTLA = [1 1 1];
    else
        MatOut = scaleMat;    scaleOut = scaleCond;
    end;

    [~, s, w] = svd(MatOut,0);
    
    if svdFlag(3) == 1
        permuteOrder = 1:size(Mat,2);
        fig = fig + 1;
        setFigure(fig, 2, figPosition, {figNamePrefix; 'SVD after scaling'}, ...
                  figXlabel, figYlabel, [], [], []);
        plotSVD(fig, diag(s), w, figRHS, permuteOrder, plotEigVec, [], [], [], sizeSV, widthSV);
        saveFigure(fig, fullfile('Figures', [figNamePrefix, 'svd-03']), saveFigFlag);
    end;
        
    if svdFlag(4) == 1
        permuteOrder = enhanceDiagDom(w);
        [~, s, w] = svd(MatOut(:,permuteOrder),0);
        fig = fig + 1;
        setFigure(fig, 2, figPosition, {figNamePrefix; 'Sorted SVD after scaling'}, ...
                  figXlabel, figYlabel, [], [], []);
        plotSVD(fig, diag(s), w, figRHS, permuteOrder, plotEigVec, [], [], [], sizeSV, widthSV);
        saveFigure(fig, fullfile('Figures', [figNamePrefix, 'svd-04']), saveFigFlag);
    end;
    
else
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> WARNING: No optimization of the condition number \n'); 
    fprintf('    To optimize the condition number, make svdFlag(3) = 1 or svdFlag(4) = 1 \n');
    fprintf('    and rerun the script -- PAUSE \n');   pause;
end;

fignext = fig + 1;

end    % of the function