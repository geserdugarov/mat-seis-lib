%% *plotDir*
% Plot directions on the unit sphere

%%
% *Input:*

% fig            - [scalar] number of an existing figure
% vec            - [3, :] array of 3D vectors to plotted
% vecGuide       - [3, :] array of 3D vectors serving as a guide: 
%                  if dot(vec(:,i), vecGuide(:,i)) < 0, the sign of vec(:,i) is changed
% vecMarkerShape - [char] marker shape 
% vecMarkerSize  - [scalar] marker size  
% vecMarkerColor - [1, 3] color array for the markers
% lightFlag      - [scalar] equal to 1 (default) or 0 turn on or off Matlab function 'camlight'

%%
% *Output:* none

%%
% *Author:* Vladimir Grechka 2012 - 2014

%%
function plotDir(fig, vec, vecGuide, vecMarkerShape, vecMarkerSize, vecMarkerColor, lightFlag)
%% Settings, checks, and defaults
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
% narginTrue = nargin(thisFileName);   
% narginchk(narginTrue, narginTrue);

grayColor1 = [1 1 1];  grayColor2 = 0.8*grayColor1;

if nargin < 7;   lightFlag = 1;   end;

if isempty(vecGuide) == 1;  vecGuide = vec;   end;

if length(vecMarkerColor) ~= 3
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Variable ''vecMarkerColor'' is supposed to be [1, 3] vector \n');
    fprintf('>>> Instead, its value is \n');
    display(vecMarkerColor);
      error('>>> STOP');
end;

% Check the vector sizes
if size(vec,1) ~= 3  ||  size(vecGuide,1) ~= 3 || size(vec,2) ~= size(vecGuide,2)
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Sizes of arrays ''vec'' and/or ''vecGuide'' are incorrect \n');
    size(vec)
    size(vecGuide)
      error('>>> STOP');
end;

%% Plot the quadrants 
scf0 = 1.0001;  scf1 = 1.01;
iangle = u2u(0:1:360, 'deg', 'rad');  
v1 = scf0*sin(iangle);   v2 = scf0*cos(iangle);   v3 = zeros(size(v1)); 

figure(fig);
plot3(v1, v2, v3, '-', 'LineWidth', 3, 'Color', grayColor2);
plot3(v2, v3, v1, '-', 'LineWidth', 3, 'Color', grayColor2);
plot3(v3, v1, v2, '-', 'LineWidth', 3, 'Color', grayColor2);

%% Plot the unit sphere
colormap(grayColor1);
[X,Y,Z] = sphere(36);
surf(X,Y,Z, 'FaceColor', 'interp', 'EdgeColor', grayColor2, 'FaceLighting', 'phong');
 
%% Plot the directions
vecAll = zeros(3,0);   
for i = 1:size(vec,2)
    % Check the array sizes
    if sum(isnan(vec(:,i))) == 0
        vec1 = vec(:,i);
   
        if norm(vec1) > eps
            % Change the direction of vec1 if necessary
            if dot(vec1, vecGuide(:,i)) < 0;   vec1 = -vec1;   end;
            vec1 = scf1*vec1/norm(vec1);
            vecAll = cat(2, vecAll, vec1);
            
            % Plot vec1
            plot3(vec1(1), vec1(2), vec1(3), vecMarkerShape, 'LineWidth', 0.5, ...
                'MarkerSize', vecMarkerSize, ...
                'MarkerEdgeColor', 'k', 'MarkerFaceColor', vecMarkerColor);
%                'MarkerEdgeColor', vecMarkerColor, 'MarkerFaceColor', vecMarkerColor);
%                'MarkerEdgeColor', 1 - vecMarkerColor, 'MarkerFaceColor', vecMarkerColor);
        else
            fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
            fprintf('>>> Warning: Vector #%g has zero length \n', i);
        end;
   end;
end;

%% Make proper view
vecAll(any(isnan(vecAll), 2), :) = [];  % remove NaN's from 'vecAll'
vecSum = sum(vecAll, 2); 
vecSum = vecSum/norm(vecSum);
[azimuth, elevation, ~] = cart2sph(vecSum(1), vecSum(2), vecSum(3));
viewAzim = u2u(azimuth + pi/2, 'rad', 'deg');  % pi/2 is added because the frame is left-handed
viewElev = u2u(-elevation,     'rad', 'deg');  % minus sign accounts for the Z-axis pointing down
if isnan(viewAzim) == 0  &&  isnan(viewElev) == 0  && ...
   isinf(viewAzim) == 0  &&  isinf(viewElev) == 0  
    view(viewAzim, viewElev);
    if lightFlag == 1;   camlight(-viewAzim, viewElev, 'infinite');   end;
else                     % the minus sign in '-viewAzim' directs light from the side
    view(3);   
    if lightFlag == 1;   camlight('headlight', 'infinite');   end;
end;

axis('vis3d');   daspect([1 1 1]);   % axis('tight');    
drawnow;

end    % of the function
