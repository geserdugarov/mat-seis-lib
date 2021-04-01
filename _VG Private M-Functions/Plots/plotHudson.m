%% *plotHudson*
% Make Hudson plot of seismic moment tensors  

%%
% *Input:*

% fig      - [scalar] number of an existing figure
% MTI      - 'MTI' structure produced by function 'smtAttributes.m'
% plotType - [scalar] equal to 'diamond' or 'skewed' defining the plot type
% plotting - structure that controls the plot appearance

%%
% *Output:* none

%%
% *Author:* Vladimir Grechka 2014 - 2015

%% 
% *References:* 
%
% * Hudson, J. A., R. G. Pearce, R. M. Rrogers, 1989, Source type plot for inversion of the moment 
%       tensor: Journal of Geophysical Research, 94, no. B1, 765-774.
%
% * Vavrycuk, V., 2014, Moment tensor decompositions revisited: Journal of Seismology, 18, no. 4,
%       ISSN 1383-4649, DOI 10.1007/s10950-014-9463-y.

%%
function [] = plotHudson(fig, MTI, plotType)
%% Settings and defaults
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

%% Set up the figure
figure(fig);  lightGray = 0.9*[1 1 1];  
set(gca, 'YDir', 'normal');      
grid off;  axis equal;  axis off;
plot([0, 0], [-1, 1], 'k-', 'LineWidth', 2);  % ISO axis

plotting = struct;
plotting.annotationSize = 8;  
plotting.eventMarkerSize = 10;  
plotting.eventMarkerShape = 's';  
plotting.fontSize = 14;

%% Annotate Hudson plot
if strcmp(plotType, 'diamond') == 1
    % Diamond tau-k plot
    set(gca, 'XDir', 'reverse');    
    plot([-1, 0, 1, 0, -1], [0, -1, 0, 1, 0], 'k-', 'LineWidth', 2);
    plot([-1, 1], [0, 0], 'k-', 'LineWidth', 2);
    text(-1.04, 0, '-CLVD', 'HorizontalAlignment', 'left', 'VerticalAlignment',  'middle', ...
        'FontSize', plotting.fontSize, 'FontName', 'Times New Roman');
    text( 1.03, 0, '+CLVD', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
        'FontSize', plotting.fontSize, 'FontName', 'Times New Roman');
    xlim([-1.1, 1.1]);  ylim([-1.1, 1.1]);
    
    % Color legend
    text(1.1, -0.6, 'Color legend:', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
        'FontSize', plotting.fontSize, 'FontName', 'Times New Roman');
    text(1.1, -0.72, ['{\color{red}ISO, ', '\color{green}CLVD, ', '\color{blue}DC}'], ...  
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
        'FontSize', plotting.fontSize, 'FontName', 'Times New Roman');

elseif strcmp(plotType, 'skewed') == 1
    % Skewed diamond, equal-area u-nu plot
    plot([-4/3, 0, 4/3, 0, -4/3], [-1/3, -1, 1/3, 1, -1/3], 'k-', 'LineWidth', 2);
    plot([-4/3, 4/3], [-1/3, 1/3], 'k-', 'LineWidth', 2);
    text(-1.06,  0, '+CLVD', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
        'FontSize', plotting.fontSize, 'FontName', 'Times New Roman');
    text( 1.06,  0, '-CLVD', 'HorizontalAlignment', 'left',  'VerticalAlignment', 'middle', ...
        'FontSize', plotting.fontSize, 'FontName', 'Times New Roman');
    xlim([-1.5, 1.5]);  ylim([-1.1, 1.1]);

    % Color legend
    text(-1.3, -0.8, 'Color legend:', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
        'FontSize', plotting.fontSize, 'FontName', 'Times New Roman');
    text(-1.3, -0.92, ['{\color{red}ISO, ', '\color{green}CLVD, ', '\color{blue}DC}'], ...  
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
        'FontSize', plotting.fontSize, 'FontName', 'Times New Roman');
end;

%% Plot moment/source tensors on Hudson plot
if strcmp(plotType, 'diamond') == 1
    % Diamond tau-k plot
    for i = 1:size(MTI.mISO, 1)
        colorMTI = [abs(MTI.fracISO(i,1)), abs(MTI.fracCLVD(i,1)), abs(MTI.fracDC(i,1))];
        plot(MTI.fracCLVD(i,1), MTI.fracISO(i,1), plotting.eventMarkerShape, 'LineWidth', 0.5, ...
            'MarkerSize', plotting.eventMarkerSize, ...
            'MarkerFaceColor', colorMTI, 'MarkerEdgeColor', lightGray);
    end;

elseif strcmp(plotType, 'skewed') == 1
    for i = 1:size(MTI.mISO, 1)
        % Restore the original eigenvalues of the moment tensor (Vavrucyk, 2014, eqs A1, A2)
        M(2) = MTI.M0(i,1)*(MTI.fracISO(i,1) - MTI.fracCLVD(i,1)/2);
        if MTI.fracCLVD(i,1) >= 0
            M(1) = MTI.M0(i,1)*(MTI.fracISO(i,1) + MTI.fracDC(i,1) + MTI.fracCLVD(i,1));
            M(3) = MTI.M0(i,1)*(MTI.fracISO(i,1) - MTI.fracDC(i,1) - MTI.fracCLVD(i,1)/2);
        else
            M(1) = MTI.M0(i,1)*(MTI.fracISO(i,1) + MTI.fracDC(i,1) - MTI.fracCLVD(i,1)/2);
            M(3) = MTI.M0(i,1)*(MTI.fracISO(i,1) - MTI.fracDC(i,1) + MTI.fracCLVD(i,1));
        end;
        % Renormalize the eigenvalues (Vavrucyk, 2014, eqs A1, A2)
        mM = M/max(abs(M));
        u  = (-2/3)*(mM(1) + mM(3) - 2*mM(2));
        v  =  (1/3)*(mM(1) + mM(2) + mM(3));
        
        colorMTI = [abs(MTI.fracISO(i,1)), abs(MTI.fracCLVD(i,1)), abs(MTI.fracDC(i,1))];
        plot(u, v, plotting.eventMarkerShape, 'LineWidth', 0.5, ...
            'MarkerSize', plotting.eventMarkerSize, ...
            'MarkerFaceColor', colorMTI, 'MarkerEdgeColor', lightGray);
    end;

else
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Incorrect value of variable plotType = %s \n', plotType);
    fprintf('    Its legitimate values are ''diamond'' or ''skewed'' -- PAUSE \n');  pause;  
end;

%% More plot annotation
if strcmp(plotType, 'skewed') == 1
    % Tensile crack points
    plot(-4/9,  5/9, 'o', 'LineWidth', 0.5, 'MarkerSize', plotting.annotationSize, ...
        'MarkerFaceColor', lightGray, 'MarkerEdgeColor', 'k');
    text(-4/9-0.05,  5/9 + 0.05, 'Tensile crack', 'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'middle', 'FontSize', plotting.fontSize, 'FontName', 'Times New Roman');
    text(-4/9-0.25,  5/9 - 0.05, 'opening', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', 'FontSize', plotting.fontSize, 'FontName', 'Times New Roman');

    plot( 4/9, -5/9, 'o', 'LineWidth', 0.5, 'MarkerSize', plotting.annotationSize, ...
        'MarkerFaceColor', lightGray, 'MarkerEdgeColor', 'k');
    text( 4/9+0.11, -5/9 + 0.05, 'Tensile crack', 'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'middle', 'FontSize', plotting.fontSize, 'FontName', 'Times New Roman');
    text( 4/9+0.30, -5/9 - 0.05, 'closure', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', 'FontSize', plotting.fontSize, 'FontName', 'Times New Roman');
end;

% ISO points
plot(0,  1, 'o', 'LineWidth', 0.5, 'MarkerSize', plotting.annotationSize, ...
    'MarkerFaceColor', lightGray, 'MarkerEdgeColor', 'k');
text(0,  1.02, 'Explosion', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
    'FontSize', plotting.fontSize, 'FontName', 'Times New Roman');
plot(0, -1, 'o', 'LineWidth', 0.5, 'MarkerSize', plotting.annotationSize, ...
    'MarkerFaceColor', lightGray, 'MarkerEdgeColor', 'k');
text(0, -1.02, 'Implosion', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
    'FontSize', plotting.fontSize, 'FontName', 'Times New Roman');

% CLVD points
plot(-1,  0, 'o', 'LineWidth', 0.5, 'MarkerSize', plotting.annotationSize, ...
    'MarkerFaceColor', lightGray, 'MarkerEdgeColor', 'k');
plot( 1,  0, 'o', 'LineWidth', 0.5, 'MarkerSize', plotting.annotationSize, ...
    'MarkerFaceColor', lightGray, 'MarkerEdgeColor', 'k');

drawnow;
   
end    % of the function

%%