%% *plotLcurve*
% Compute and plot the L-curve using the 'Regularization Tools' package developed by P.C.Hansen 
% (<http://www.mathworks.com/matlabcentral/fileexchange/52>)

%%
% *Input:*

% fig           - [scalar] number of the figure to be created
% errData       - [scalar] relative data error
% F             - [:, :] arrays of scaled Frechet derivatives
% model         - [1, size(F,2)] model-parameter vector
% figPosition   - [4, 1] vector specified the figure position and size
% figNamePrefix - [char] prefix of the file name for saving generated figures
% markerSize    - [scalar] marker size for plotting the L-curve
% lineWidth     - [scalar] line width for plotting the L-curve
% saveFigFlag   - [3, 1] array containing ones and zeros used in 'saveFigure'

%%
% *Output:* none

%%
% *Author:* Vladimir Grechka 2012 2013

%%
function plotLcurve(fig, errData, F, model, ...
    figPosition, figNamePrefix, markerSize, lineWidth, saveFigFlag)             
%% Settings and defaults
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

global letterDrive

if isempty(markerSize) == 1;    markerSize = 4;      end; 
if isempty(lineWidth)  == 1;    lineWidth  = 0.5;    end; 

%% Compute the data perturbation due to 'timeTrue' and add noise
data    = F*model;
dataErr = data.*(1 + errData*randn(size(data)));

%% Apply Hansen's package to compute the L-curve
folderNameHansen = [letterDrive, ...
    '\Users\erf\Grechka\mSeismic Tools\_Private M-Functions\Regularization Tools by PCH\regu'];
addpath(folderNameHansen);
[U, s, ~] = csvd(F);
[regCorner, rho, eta, regParam] = l_curve_VG(U, s, dataErr, 'Tikh');
rmpath(folderNameHansen);

iCorner = find(abs(regParam - regCorner) == min(abs(regParam - regCorner)));

%% Plot the L-curve
figure(fig);    % 'setFigure' cannot be used because it makes linear rather than log-log scale
set(fig, 'Color', [1 1 1]);     % set color of the figure window
set(fig, 'Position', figPosition);    
box on;                         % plot border around the figure
hold off;  
loglog(rho, eta, 'o', 'LineWidth', lineWidth, 'MarkerSize', markerSize, ...
    'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
hold on;    grid on;

text(rho(iCorner), eta(iCorner), [' L-corner = ' num2str(regCorner)], ...
    'VerticalAlignment', 'bottom', 'Fontsize', 14, 'Fontname', 'Times', ...
    'BackgroundColor', 'none');
plot(rho(iCorner), eta(iCorner), '^', 'LineWidth', lineWidth, 'MarkerSize', markerSize + 7, ...
    'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'w');  axis('tight');
set(gca, 'Fontsize', 14, 'Fontname', 'Times'); 
xlabel('Log_{10} of residual norm  || {\itF} {\itm}_{\lambda} - {\itd}_{ }||');   
ylabel('Log_{10 } of solution norm  ||  {\it{m}}_{\lambda }||');   
title({figNamePrefix; 'L-curve'}, 'Fontsize', 16, 'Fontname', 'Times', 'Interpreter', 'none');
drawnow;

saveFigure(fig, fullfile('Figures', [figNamePrefix, 'L-curve']), saveFigFlag);

%% Result of the L-curve analysis
fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
fprintf('>>> Tikhonov''s L-corner (= %g) is located \n', regCorner);  
if iCorner == length(regParam);  
    fprintf('    at the end of the regularization interval -> no regularization needed \n \n'); 
else
    fprintf('    in the middle of the regularization interval -> regularization could be helpful... \n');  
%    fprintf('>>> Proceeding without reqularization... \n \n'); 
end;

end    % of the function