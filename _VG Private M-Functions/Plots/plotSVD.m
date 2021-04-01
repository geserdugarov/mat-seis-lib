%% *plotSVD*
% Plot singular values and eigenvectors resulting from the SVD  

%%
% *Input:*

% fig             - [scalar] number of the existing figure   
% s               - array of the singular values
% w               - the eigenvector matrix
% varCaption      - [cell] array of captions 
% permutation     - [1, :] permutation order of array varCaption 
% plotEigVec      - [scalar] parameter that controls plotting the eigenvalue matrix
%                   (*) Both the eigenvalue matrix and the singular values are plotted if 
%                       length(s) <= plotEigVec; only singular values are plotted otherwise
% flagCaption     - [scalar] flag equal to 1 or 0 to indicate whether to display the captions
% fontSizeCaption - [scalar] font size for the captions
% dxticks         - [scalar] increment of labeling the x-axis
% sizeSV          - [scalar] size of circles indicating the singular values
% widthSV         - [scalar] width of line connecting the singular values

%%
% *Output:* none

%%
% *Author:* Vladimir Grechka 2012 - 2014

%%
function plotSVD(fig, s, w, varCaption, permutation, plotEigVec, ...
                 flagCaption, fontSizeCaption, dxticks, sizeSV, widthSV)
%% Settings and defaults
[~, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

tolSVD = 1.e-16;    
lightgray = 0.3*[1 1 1]; 

if isempty(plotEigVec) == 1;       plotEigVec = 50;        end; 
if isempty(flagCaption) == 1;      flagCaption = 1;        end; 
if isempty(fontSizeCaption) == 1;  fontSizeCaption = 16;   end; 
if isempty(dxticks) == 1
    dxticks0 = 5:5:1000;  
    dxticks  = dxticks0(abs(dxticks0 - floor(length(s)/6)) == ...
                    min(abs(dxticks0 - floor(length(s)/6))) );    
end; 
if isempty(sizeSV)  == 1;    sizeSV = 6;    end; 
if isempty(widthSV) == 1;   widthSV = 2;    end; 

%% Normalize the singular values  
sNorm  = s/max(s);   
sNorm(sNorm < tolSVD) = tolSVD;      % replace small singular values with tolSVD
noSV  = length(sNorm);   noSV2    = noSV^2;
logSV = log10(sNorm);    minLogSV = min(logSV);

%% Setup appearance of the figure
figure(fig);
%set(fig, 'Color', 'w');             % set the background color to white

% Map the eigenvector interval [noSV + 1/2 : 1/2] onto the logSV interval [0 : minLogSV] 
% [1/2 : noSV + 1/2] = a*[0 : minLogSV] + b
a = noSV/minLogSV;    b = 1/2;
splot = a*logSV + b;

% Define yticks
if logSV(end) < -10.0;   
    dyticks = 2.0;
elseif logSV(end) <  -5.0;   
    dyticks = 1.0;
elseif logSV(end) <  -2.5    
    dyticks = 0.5;
else
    dyticks = 0.2;      
end;
yticks = sort(a*(dyticks*floor(min(logSV)/dyticks) : dyticks : 0) + b, 'ascend');   

% Add labels
set(gca, 'YTick', yticks); 
set(gca, 'YTickLabel', sort(dyticks*floor(min(logSV)/dyticks) : dyticks : 0, 'descend')); 

%%
if noSV > plotEigVec
    %% Plot only singular values
    plot(1:noSV, splot, '-', 'LineWidth', widthSV, 'Color', 'k');     
    plot(1:noSV, splot, 'o', 'LineWidth', 0.5, 'MarkerSize', sizeSV, ...
        'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');     
else
    %% Plot the eigenvector matrix
    set(gca, 'XTick',      (dxticks : dxticks : noSV)); 
    set(gca, 'XTickLabel', (dxticks : dxticks : noSV)); 

    % Create arrays for Matlab function 'patch'  
    x0 = [0.5, 1.5, 1.5, 0.5]';
    y0 = [0.5, 0.5, 1.5, 1.5]';
    xdata = NaN(4, noSV2);    ydata = NaN(4, noSV2);    cdata = NaN(1, noSV2); 
    for ix = 1:noSV;    
        for iy = 1:noSV
            xdata(:, noSV*(iy-1)+ix) = x0 + (ix - 1); 
            ydata(:, noSV*(iy-1)+ix) = y0 + (iy - 1); 
            cdata(1, noSV*(iy-1)+ix) = 1 - abs(w(iy,ix));
        end;   
    end;
    grid off;    hold on;
    colormap('gray'); 
    p = patch(xdata, ydata, cdata);   set(p, 'EdgeColor', 'k');
        
    %% ... and the singular values
    splot(1) = splot(1) - 1.e-3;  % adjust splot(1) to move in the frame of the plot
    plot(1:noSV, splot, '-', 'LineWidth', widthSV+2, 'Color', 'w'); 
    plot(1:noSV, splot, '-', 'LineWidth', widthSV,   'Color', lightgray); 
    plot(1:noSV, splot, 'o', 'LineWidth', 1, 'MarkerSize', sizeSV, ...
        'MarkerEdgeColor', 'w', 'MarkerFaceColor', lightgray, 'Color', lightgray); 
    axis('image');

    % Label the rows
    if flagCaption == 1
        for ij = 1:size(varCaption,2)
            xshift = 1;
            if mod(ij, 2) == 0  &&  noSV > 30
                xshift = 1 + 0.07*noSV;     % shit of text in plotted squares
            end;
            text(noSV+xshift, ij+1/4, varCaption{1,permutation(ij)}, ...
                'Fontsize', fontSizeCaption, 'Fontname', 'Times', 'VerticalAlignment', 'Middle');
        end;
    end;
end;
drawnow;

end    % of the function