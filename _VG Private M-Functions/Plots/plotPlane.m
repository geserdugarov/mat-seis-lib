%% *plotPlane*
% Plot a plane 

%%
% *Input:*

% fig        - [scalar] number of an existing figure
% corners    - [1, 6] array of coordinates of the box corners 
%              [xmin, xmax, ymin, ymax, zmin, zmax], such as in Matlab function 'axis'
%              (*) A plane will be plotted in a box specified by the corners
% plane      - [5, :] array specified by function 'setInterface' or 'setFault'
% inputUnits - [string] input length units
% plotUnits  - [string] length units for plotting 
% planeColor - [1, 3] color array
% opacity    - [scalar] parameter ranging from 0 (transparent) to 1 (opaque) to control
%              opacity through Matlab function 'alpha'

%%
% *Output:* none

%%
% *Author:* Vladimir Grechka 2012 - 2014

%%
function plotPlane(fig, corners, plane, inputUnits, plotUnits, planeColor, opacity)
%% Settings and checks
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);
    
noCell = 10;    % the number of rectangles to represent a plane

if size(plane, 2) == 0  
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Empty ''plane'' array \n');   
    return;
else
    % Array plane is not empty: make grids for plotting the planes
    if isempty(inputUnits) == 0
        plane(1:3,:) = u2u(plane(1:3,:), inputUnits, plotUnits, 0);
    end;
    dx1 = (corners(2) - corners(1))/noCell;    
    dx2 = (corners(4) - corners(3))/noCell;    
    
    % Correct for the plane dip
    sinMaxDip = max(sqrt(1 - plane(4,:).^2 - plane(5,:).^2));
    if sinMaxDip < 0.01    
        sinMaxDip = 0.01;    
    end;
    dx = sinMaxDip*max([dx1, dx2]);
    
    % Avoid zero dx
    if dx < 1.e-6;  
        dx = 0.5*(max(plane(3,:)) - min(plane(3,:)));   
    end;
    x1 = corners(1) - dx : dx : corners(2) + dx;
    x2 = corners(3) - dx : dx : corners(4) + dx;

    noPln = size(plane,2);
    x3 = NaN(length(x2), length(x1), noPln);
    for ix1 = 1:length(x1)
        for ix2 = 1:length(x2)
            x3(ix2,ix1,:) = planeDepth(plane, [x1(ix1), x2(ix2)]', 1:noPln);
        end;
    end;
end;

% Check whether the interface depths are in the box
if max(max(max(x3))) < corners(5) || min(min(min(x3))) > corners(6)
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> WARNING: The planes are plotted are outside the specified box \n \n');    
%    return;
end;
 
%% Plotting
figure(fig);  
axis(corners);    grid on;   hold on;

lightColor1 = 0.5 + 0.5*planeColor;    
lightColor2 = 0.9 + 0.1*planeColor;
colormap([lightColor1; lightColor2]);

for ipln = 1:noPln
    if isinf(plane(3,ipln)) == 0
        surf(x1, x2, x3(:,:,ipln), 'LineWidth', 0.5, ...
            'EdgeColor', lightColor1, 'FaceColor', lightColor2);   hold on;  
    end;
end;
alpha(opacity);

axis('image');   axis('tight');   
view(3);   drawnow;

end    % of the function