%% *plotBB*
% Plot seismic moment tensors as beach ball diagrams

%%
% *Input:*

% fig       - [scalar] number of an existing figure
% angles    - [:, 3] array of the fault strikes, dips, and rakes (all in degrees) 
% xSou      - [dim, :] array of the event locations
%             (*) dim = 2 - 2D plot
%             (*) dim = 3 - 3D plot
% sizeBB    - [1, size(xSou,2)] array of radii of the plotted beach balls
% colorPres - [:, 3] color array for the pressure portion of the beach balls
% colorTens - [:, 3] color array for the tension portion of the beach balls
% alphaNum  - [scalar] number determining transparancy in Matlab function 'alpha'

%%
% *Output:* 

% nullAxis  - [2, :] array of coordinates of the null axes defined as the intersections of 
%             the primary and auxiliary focal planes

%%
% *Author:* Vladimir Grechka 2012 - 2014

%%
function [nullAxis] = plotBB(fig, angles, xSou, sizeBB, colorPres, colorTens, alphaNum)
%% Settings and defaults
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

dim = size(xSou, 1);
dphi = 3;   phi = (0 : dphi : 360);     

if size(angles, 1) ~= size(xSou, 2)
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Inconsistent numbers of seismic moment tensors (= %g) and \n', size(angles,1));
    fprintf('    event locations (%g) \n \n', size(xSou,2));    
      error('>>> STOP');
end;

figure(fig);

%% Loop over events
nullAxis = NaN(2, size(angles,1));
for ievnt=1:size(angles,1)

    perim = [sind(phi); cosd(phi)];                     % perimeter of the unit cirle

    if angles(ievnt,2) < 0  ||  angles(ievnt,2) > 90
        fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
        fprintf('>>> Dip = %g (deg) for event %g is outside the plausible range [0, 90] \n', ...
                angles(ievnt,2), ievnt);    
        fprintf('>>> The event is skipped -- PAUSE \n');    pause;

    else
    
        if abs(angles(ievnt,2)) == 90 
            fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
            fprintf('>>> The dip of %g degrees for event %g entails ambiguity in coloring \n', ...
                    angles(ievnt,2), ievnt);    
            fprintf('    its beach ball because the foot-wall and the hanging-wall of a fault \n');
            fprintf('    are indeterminable -- PAUSE \n');    pause;           
%            angles(ievnt,2) = angles(ievnt,2) + 1; % ... perhaps correct the angle
        end;

        if size(sizeBB, 2) < ievnt
            fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
            fprintf(['>>> WARNING: Undefined size of beach ball for event %g is replaced', ...
                     ' with that for event 1 \n'], ievnt);
            sizeBB(1,ievnt) = sizeBB(1,1);
        end;

        % Place the rake in the range [-180, 180] 
        angles(ievnt,3) = mod(angles(ievnt,3), 360);
        if angles(ievnt,3) > 180;   angles(ievnt,3) = angles(ievnt,3) - 360;   end;
    
        % Get the nodal planes
        [plane1, plane2] = smtNL(angles(ievnt, 1:3));
        
        %------------------------------------------------------------------------------------------
        % Find intersections of the nodal planes with the outer circle
        q1 = u2u([atan2(plane1(1,1), plane1(2,1)), atan2(plane1(1,end), plane1(2,end))], ...
                    'rad', 'deg');
        q2 = u2u([atan2(plane2(1,1), plane2(2,1)), atan2(plane2(1,end), plane2(2,end))], ...
                    'rad', 'deg');
        for i=1:2
            if q1(i) < 0;   q1(i) = q1(i) + 360;    end;
            if q2(i) < 0;   q2(i) = q2(i) + 360;    end;
        end;
    
        % Remove columns of NaNs
        plane1(:,any(isnan(plane1),1)) = []; 
        plane2(:,any(isnan(plane2),1)) = [];    

        % Separate the focal sphere into quadrants or hemispheres
        if isempty(plane1) == 1    % the fault plane is not visible -- two hemispheres
            if q2(2) > q2(1)
                perimArray1 = [(q2(2) :  dphi : 360 - dphi/2), (0 : dphi : q2(1) - dphi/2), q2(1)];
                perimArray2 = [(q2(2) : -dphi : q2(1) + dphi/2), q2(1)];
            else
                perimArray1 = [(q2(2) :  dphi : q2(1) - dphi/2), q2(1)];
                perimArray2 = [(q2(2) : -dphi : dphi/2), (360 : -dphi : q2(1) + dphi/2), q2(1)];
            end;
    
            % Ad-hoc assignment to data1 and data2 (and the corresponding coloring) because 
            % the node plane goes through the coordinate origin
            xdata1 = [plane2(1,:), sind(perimArray1)];
            ydata1 = [plane2(2,:), cosd(perimArray1)];
            cdata1 = ones(size(xdata1));
        
            xdata2 = [plane2(1,:), sind(perimArray2)];
            ydata2 = [plane2(2,:), cosd(perimArray2)];
            cdata2 = zeros(size(xdata2));

        elseif isempty(plane2) == 1    % the auxiliary plane is not visible -- two hemispheres
            if q1(2) > q1(1)
                perimArray1 = [(q1(2) :  dphi : 360 - dphi/2), (0 : dphi : q1(1) - dphi/2), q1(1)];
                perimArray2 = [(q1(2) : -dphi : q1(1) + dphi/2), q1(1)];
            else
                perimArray1 = [(q1(2) :  dphi : q1(1) - dphi/2), q1(1)];
                perimArray2 = [(q1(2) : -dphi : dphi/2), (360 : -dphi : q1(1) + dphi/2), q1(1)];
            end;

            % Ad-hoc assignment to data1 and data2 (and the corresponding coloring) because 
            % the node plane goes through the coordinate origin
            xdata1 = [plane1(1,:), sind(perimArray1)];
            ydata1 = [plane1(2,:), cosd(perimArray1)];
            cdata1 = ones(size(xdata1));
            
            xdata2 = [plane1(1,:), sind(perimArray2)];
            ydata2 = [plane1(2,:), cosd(perimArray2)];
            cdata2 = zeros(size(xdata2));

        else
            ind1 = 0;   ind2 = 0;   dist = 2;   % both the fault and auxiliary planes are visible
            for i1=1:size(plane1,2)-1           % find the null axis (N-axis)
                tmp1 = (plane1(:,i1+1) + plane1(:,i1))/2; 
                for i2=1:size(plane2,2)-1
                    tmp2 = (plane2(:,i2+1) + plane2(:,i2))/2; 
                    dd = sqrt( dot(tmp2 - tmp1, tmp2 - tmp1) );
                    if dd < dist
                        dist = dd;
                        ind1 = i1;   ind2 = i2;
                    end;
                end;
            end;
            
            k1 = (plane1(2,ind1+1) - plane1(2,ind1))/(plane1(1,ind1+1) - plane1(1,ind1));
            k2 = (plane2(2,ind2+1) - plane2(2,ind2))/(plane2(1,ind2+1) - plane2(1,ind2));
            b1 = plane1(2,ind1) - k1*plane1(1,ind1);
            b2 = plane2(2,ind2) - k2*plane2(1,ind2);
            axisN(1) = (b1 - b2)/(k2 - k1);         % N-axis as solution of linear equations
            axisN(2) = (k2*b1 - k1*b2)/(k2 - k1);   % (Korn and Korn, 1984, eq 2.3-4)     
            nullAxis(:,ievnt) = [axisN(1); axisN(2)];
            
            % Make the quadrants
            % Quadrant 1: axisN - q1(1) - q2(1) - axisN
            clear perimArray
            perimArray = smtPerimAzim(q1(1), q2(1), dphi);
            xdataQ1 = [axisN(1), plane1(1, (ind1 : -1 : 1)), sind(perimArray), ...
                                 plane2(1, (1 : ind2)), axisN(1)];
            ydataQ1 = [axisN(2), plane1(2, (ind1 : -1 : 1)), cosd(perimArray), ...
                                 plane2(2, (1 : ind2)), axisN(2)];
            cdataQ1 = ones(size(xdataQ1));

            % Quadrant 2: axisN - q2(1) - q1(2) - axisN
            clear perimArray
            perimArray = smtPerimAzim(q2(1), q1(2), dphi);
            xdataQ2 = [axisN(1), plane2(1, (ind2 : -1 : 1)), sind(perimArray), ...
                                 plane1(1, (end  : -1 : ind1 + 1)), axisN(1)];
            ydataQ2 = [axisN(2), plane2(2, (ind2 : -1 : 1)), cosd(perimArray), ...
                                 plane1(2, (end  : -1 : ind1 + 1)), axisN(2)];
            cdataQ2 = ones(size(xdataQ2));
        
            % Quadrant 3: axisN - q1(2) - q2(2) - axisN
            clear perimArray
            perimArray = smtPerimAzim(q1(2), q2(2), dphi);
            xdataQ3 = [axisN(1), plane1(1, (ind1 + 1 : end)), sind(perimArray), ...
                                 plane2(1, (end : -1 : ind2 + 1)), axisN(1)];
            ydataQ3 = [axisN(2), plane1(2, (ind1 + 1 : end)), cosd(perimArray), ...
                                 plane2(2, (end : -1 : ind2 + 1)), axisN(2)];
            cdataQ3 = ones(size(xdataQ3));
            
            % Quadrant 4: axisN - q1(1) - q2(2) - axisN
            clear perimArray
            perimArray = smtPerimAzim(q1(1), q2(2), dphi);
            xdataQ4 = [axisN(1), plane1(1, (ind1 : -1 : 1)), sind(perimArray), ...
                                 plane2(1, (end  : -1 : ind2 + 1)), axisN(1)];
            ydataQ4 = [axisN(2), plane1(2, (ind1 : -1 : 1)), cosd(perimArray), ...
                                 plane2(2, (end  : -1 : ind2 + 1)), axisN(2)];
            cdataQ4 = ones(size(xdataQ4));
        
            % Color the quandrants
            % Find the quadrant containing the coordinate origin
            gravityCenter = [[mean(xdataQ1), mean(xdataQ2), mean(xdataQ3), mean(xdataQ4)]; ...
                             [mean(ydataQ1), mean(ydataQ2), mean(ydataQ3), mean(ydataQ4)]];
            iQuad = find( dot(gravityCenter, gravityCenter, 1) == ...
                      min(dot(gravityCenter, gravityCenter, 1)) );
            
            % Determine the quadrant containing the tension axis
            if angles(ievnt,3) < 0;   iQuad = iQuad - 1;   end;
            if iQuad == 1  ||  iQuad == 3
                xdata1 = [xdataQ1, xdataQ3];    xdata2 = [xdataQ2, xdataQ4];    % data1 is tension
                ydata1 = [ydataQ1, ydataQ3];    ydata2 = [ydataQ2, ydataQ4];    % data2 is pressure
                cdata1 = [cdataQ1, cdataQ3];    cdata2 = [cdataQ2, cdataQ4];
            else
                xdata1 = [xdataQ2, xdataQ4];    xdata2 = [xdataQ1, xdataQ3];
                ydata1 = [ydataQ2, ydataQ4];    ydata2 = [ydataQ1, ydataQ3];
                cdata1 = [cdataQ2, cdataQ4];    cdata2 = [cdataQ1, cdataQ3];
            end;
        end;

        % Coordinate transform
        xdata1 = xSou(1,ievnt) + sizeBB(1,ievnt)*xdata1;   
        xdata2 = xSou(1,ievnt) + sizeBB(1,ievnt)*xdata2;
        ydata1 = xSou(2,ievnt) + sizeBB(1,ievnt)*ydata1;   
        ydata2 = xSou(2,ievnt) + sizeBB(1,ievnt)*ydata2;
        plane1 = repmat(xSou(1:2,ievnt), 1, size(plane1,2)) + sizeBB(1,ievnt)*plane1;   
        plane2 = repmat(xSou(1:2,ievnt), 1, size(plane2,2)) + sizeBB(1,ievnt)*plane2;
        perim  = repmat(xSou(1:2,ievnt), 1, size(perim, 2)) + sizeBB(1,ievnt)*perim;
        
        if dim == 2
            patch(xdata1, ydata1, cdata1, ...
                    'FaceColor', colorTens(ievnt,:), 'EdgeColor', colorTens(ievnt,:), ...
                    'FaceAlpha', alphaNum);
            patch(xdata2, ydata2, cdata2, ...
                    'FaceColor', colorPres(ievnt,:), 'EdgeColor', colorPres(ievnt,:), ...
                    'FaceAlpha', alphaNum);
            plot(plane1(1,:), plane1(2,:), '-', 'LineWidth', 0.5, 'Color', colorTens(ievnt,:));
            plot(plane2(1,:), plane2(2,:), '-', 'LineWidth', 0.5, 'Color', colorTens(ievnt,:));
            plot(perim(1,:),  perim(2,:),  '-', 'LineWidth', 0.5, 'Color', colorTens(ievnt,:));
        else
            zdata1 = repmat(xSou(3,ievnt), size(xdata1));
            zdata2 = repmat(xSou(3,ievnt), size(xdata2));
            patch(xdata1, ydata1, zdata1, cdata1, ...
                    'FaceColor', colorTens(ievnt,:), 'EdgeColor', colorTens(ievnt,:), ...
                    'FaceAlpha', alphaNum);
            patch(xdata2, ydata2, zdata2, cdata2, ...
                    'FaceColor', colorPres(ievnt,:), 'EdgeColor', colorPres(ievnt,:), ...
                    'FaceAlpha', alphaNum);
            plot3(plane1(1,:), plane1(2,:), xSou(3,ievnt)*ones(size(plane1(1,:))), ...
                    '-', 'LineWidth', 0.5, 'Color', colorTens(ievnt,:));
            plot3(plane2(1,:), plane2(2,:), xSou(3,ievnt)*ones(size(plane2(1,:))), ...
                    '-', 'LineWidth', 0.5, 'Color', colorTens(ievnt,:));
            plot3(perim(1,:),  perim(2,:),  xSou(3,ievnt)*ones(size( perim(1,:))), ...
                    '-', 'LineWidth', 0.5, 'Color', colorTens(ievnt,:));
        end;
%        axis('equal');   
        drawnow;
    end;    % of the first 'if' statement in the loop

end;    % of the 'ievnt' loop
    
end    % of the function