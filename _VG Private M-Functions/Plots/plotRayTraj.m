%% *plotRayTraj*
% Plot ray trajectories

%%
% *Input:*

% fig        - [scalar] number of an existing figure
% ray        - [3, :, :] array of the ray trajectories produced by 'art3D'
% inputUnits - [string] input length units
% plotUnits  - [string] length units for plotting 
% noSou      - [scalar] number of sources
% noRec      - [scalar] number of receivers
% noWave     - [scalar] number of wave codes

%%
% *Output:* none

%%
% *Author:* Vladimir Grechka 2012 

%%
function plotRayTraj(fig, ray, inputUnits, plotUnits, noSou, noRec, noWave)
%% Settings  
[~, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

rayColor  = setPlottingColors;
noColor = size(rayColor,1);

%% Plotting
figure(fig);
ray = u2u(ray, inputUnits, plotUnits, 0);           % change the units 
for iwave = 1:noWave
    colorIndex = mod(iwave, noColor);               % cycle the ray colors
    if colorIndex == 0;   colorIndex = noColor;   end;
    
    for isou = 1:noSou
        for irec = 1:noRec
            itime = noRec*noWave*(isou - 1) + noWave*(irec - 1) + iwave;
            plot3(ray(1,:,itime), ray(2,:,itime), ray(3,:,itime), '-', 'LineWidth', 1, ...
                'Color', rayColor(colorIndex,:));
        end;
    end;
end;

axis('image');   axis('tight');   view(3);   drawnow;

end    % of the function