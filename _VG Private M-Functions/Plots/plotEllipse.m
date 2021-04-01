%% *plotEllipse* 
% Plot 2D ellipse of 3D ellipsoid 

%%
% *Input:*

% fig      - [scalar] number of an existing figure
% x0       - [N, 1] coordinates of the center of an ellipse (N = 2) or an ellipsoid (N = 3) 
% A        - [N, N] symmetric positive definite matrix describing the ellipse in the form
%            (x - x0)'*A*(x - x0) = 1
% flag     - [scalar] equal to 1 or 0 to indicate whether axes of an ellipse are to be displayed
% colorEll - [1, 3] color array for plotting the ellipse
% widthEll - [scalar] width of line to plot the ellipse

%%
% *Output:* none

%%                    
% *Author:* Vladimir Grechka 2012 - 2014

%%
function [] = plotEllipse(fig, x0, A, flag, colorEll, widthEll)
%% Settings and defaults
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue); 
N = size(x0,1);  

if isempty(colorEll) == 1;   colorEll = 'r';   end;
if isempty(widthEll) == 1;   widthEll = 2;     end;

%% Arrays for plotting
dang = u2u(5, 'deg', 'rad');
angleAzi = (0:dang:2*pi)';   sazi = sin(angleAzi);   cazi = cos(angleAzi);
anglePol = (0:dang:  pi)';   spol = sin(anglePol);   cpol = cos(anglePol);
nAzi = length(angleAzi);   
nPol = length(anglePol);   
    
%% Translate the eigen properties of A into the semi-axes of the ellipsoid
[vv, d] = eig(A);
dd = 1./sqrt(diag(d));                  % semi-axes
dv = vv*diag(dd);                       % A = [dd(1)*vv(:,1), dd(2)*vv(:,2), dd(3)*vv(:,3)] 
    
%%
figure(fig);
if N == 2
    % Ellipse
    x = NaN(2, nAzi);
    count = 0;
    for ia=1:nAzi
        count = count + 1;
        n2 = [cazi(ia), sazi(ia)]';     % unit vector in 2D 
        x(:,count) = x0 + dv*n2;        % coordinates of a point at the ellipse
    end;  
    plot(x(1,:), x(2,:), 'o-', 'LineWidth', 0.5, ...
            'MarkerSize', 3, 'MarkerEdgeColor', colorEll, 'MarkerFaceColor', colorEll, 'Color', colorEll);
      
    if flag == 1  
        % Plot axes of the ellipse
        y1 = repmat(x0, 1, N) - dv;
        y2 = repmat(x0, 1, N) + dv;
        for i=1:N
            plot([y1(1,i), y2(1,i)], [y1(2,i), y2(2,i)], '-', 'LineWidth', widthEll, 'Color', colorEll);
        end;
    end;

elseif N == 3;
    % Ellipsoid
    y = NaN(3, nAzi, nPol);
    for ip=1:nPol;  
        for ia=1:nAzi
            n3 = [spol(ip)*cazi(ia), spol(ip)*sazi(ia), cpol(ip)]';     % unit vector in 3D 
            y(:,ia,ip) = x0 + dv*n3;        % points at the surface of ellipsoid
        end;    
    end;

    % Create arrays for 'patch'
    xdata = NaN(4, (nPol-1)*(nAzi-1));    ydata = NaN(4, (nPol-1)*(nAzi-1));
    zdata = NaN(4, (nPol-1)*(nAzi-1));    cdata = ones(1, (nPol-1)*(nAzi-1));
    count = 0;
    for ip=1:nPol-1;  
        for ia=1:nAzi-1
            count = count + 1;
            xdata(:,count) = [y(1,ia,ip), y(1,ia+1,ip), y(1,ia+1,ip+1), y(1,ia,ip+1)]';
            ydata(:,count) = [y(2,ia,ip), y(2,ia+1,ip), y(2,ia+1,ip+1), y(2,ia,ip+1)]';
            zdata(:,count) = [y(3,ia,ip), y(3,ia+1,ip), y(3,ia+1,ip+1), y(3,ia,ip+1)]';
            cdata(1,count) = 1;
        end;    
    end;
    
    % Plot transparent ellipsoid that has color specified by variable colorEll
    colormap([0 0 0]);  
    patchLine = patch(xdata, ydata, zdata, cdata, 'FaceLighting', 'phong', ...
         'FaceColor', 'interp', 'FaceAlpha', 0.05, 'EdgeColor', 'none');
    set(patchLine,'FaceColor', colorEll)

    % Remove lines that belong to the ellipsoid from the legend
    set(get(get(patchLine, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');

    % Plot the equator and main meridians
    r12 = NaN(3,nAzi);    r13 = NaN(3,nAzi);    r23 = NaN(3,nAzi); 
    for ia=1:nAzi
        n2 = [cazi(ia), sazi(ia)]';
        r12(:,ia) = x0 + [dv(:,1), dv(:,2)]*n2;
        r13(:,ia) = x0 + [dv(:,1), dv(:,3)]*n2;
        r23(:,ia) = x0 + [dv(:,2), dv(:,3)]*n2;
    end;
    plotLine = plot3(r12(1,:), r12(2,:), r12(3,:), 'Color', colorEll, 'LineWidth', widthEll);
    set(get(get(plotLine, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
    plotLine = plot3(r13(1,:), r13(2,:), r13(3,:), 'Color', colorEll, 'LineWidth', widthEll);
    set(get(get(plotLine, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
    plotLine = plot3(r23(1,:), r23(2,:), r23(3,:), 'Color', colorEll, 'LineWidth', widthEll);
    set(get(get(plotLine, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
    
    if flag == 1  
        y1 = repmat(x0, 1, N) - dv;     % construct and plot axes of an ellipsoid
        y2 = repmat(x0, 1, N) + dv;
        for i=1:N
            plot3([y1(1,i), y2(1,i)], [y1(2,i), y2(2,i)], ...
                  [y1(3,i), y2(3,i)], '-', 'LineWidth', widthEll, 'Color', colorEll);
        end;
    end;
    
else
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Incorrect ellipsoid dimension N = %g \n', N);
    fprintf('>>> The values of N should be N = 2 (ellipse) or N = 3 (ellipsoid) \n \n');
      error('>>> STOP');
end;
drawnow;

end    % of the function 
