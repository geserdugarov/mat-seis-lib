%% *plotHodogram*
% Plot the hodogram of particle motion

%%
% *Input:*

% fig              - [scalar] number of a new figure
% traceData        - [:, 3] matrix containing a single three-component trace
% flagNorm         - [scalar] normalization flag
%                    flagNorm = 0 - no normalization (default)
%                    flagNorm = 1 - normalize to the maximim of the entire seismogram
% markerSizeSample - [scalar] marker size for time samples
% markerSizeBeg    - [scalar] marker size for the beginning of the hodogram
% symbolBeg        - [char] symbol marking the beginning of the hodogram
% markerSizeEnd    - [scalar] marker size for the end of the hodogram
% symbolEnd        - [char] symbol marking the end of the hodogram
% timeBeg          - [scalar] time of the beginning of the hodogram
% timeEnd          - [scalar] time of the end of the hodogram
% timeUnits        - [string] units of time of the hodogram
% timeSample       - [scalar] time sample in 'timeUnits'
% interpFlag       - [scalar] flag governing interpolation: 0 - no, 1 - yes
% subsampleRate    - [scalar] number of samples per real sample; default(subsampleRate) = 5 

%%
% *Output:* none

%%
% *Author:* Vladimir Grechka 2012

%%
function plotHodogram(fig, traceData, flagNorm, markerSizeSample, ...
    markerSizeBeg, symbolBeg, markerSizeEnd, symbolEnd, timeBeg, timeEnd, timeUnits, timeSample, ...
    interpFlag, subsampleRate)
%% Settings and defaults
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

if isempty(traceData) == 1
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Empty trace array \n');    
    fprintf('>>> PAUSE -- Continue? \n');    
    return;
end;

if isempty(markerSizeSample) == 1;  markerSizeSample = 5;  end;   % default for 'markerSizeSample'
if isempty(markerSizeBeg) == 1;     markerSizeBeg = 24;    end;   % default for 'markerSizeBeg'
if isempty(symbolBeg) == 1;         symbolBeg = 'p';       end;   % default for 'symbolBeg'
if isempty(markerSizeEnd) == 1;     markerSizeEnd = 14;    end;   % default for 'markerSizeEnd'
if isempty(symbolEnd) == 1;         symbolEnd = 's';       end;   % default for 'symbolEnd'
if isempty(timeBeg) == 1;           timeBeg = 1;           end;   % default for 'timeBeg'
if isempty(timeEnd) == 1;           timeEnd = size(traceData, 1); % default for 'timeEnd'          
end;   
if isempty(timeUnits) == 1;         timeUnits = 'samples'; end;   % default for 'timeUnits'
if isempty(interpFlag) == 1;        interpFlag = 0;        end;   % default for 'interpFlag'

if interpFlag == 1 && isempty(subsampleRate) == 1
    subsampleRate = 5;
else
    subsampleRate = abs(subsampleRate);
end;

if subsampleRate < 1
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Variable subsampleRate < 1, it will be defaulted to 1 \n');    
    fprintf('>>> PAUSE -- Continue? \n');   pause; 
    subsampleRate = 1;
end;

%% Normalize the hodogram
if isempty(flagNorm) == 1;   flagNorm = 0;   end;                 % default for 'symbolEnd'

maxTrace = max(max(abs(traceData)));
if flagNorm == 1
    traceData = traceData/maxTrace;                               % normalize the  hodogram
    maxTrace = 1;
end;

%% Set up plotting parameters
figure(fig);
%map = flipud(colormap('jet'));
map = (colormap('jet'));
noTime = size(traceData, 1);
dcolor = (size(map, 1) - 1)/(timeSample*(noTime - 1));            % step in the color map
lightGray = 0.4*[1 1 1];

%% Plot the hodogram
if interpFlag == 1
    % Interpolate the hodogram
    tt1 = 1:size(traceData, 1);  
    tti = 1 : 1/subsampleRate : size(traceData, 1); 
    traceDataInt = zeros(length(tti), 3);
    for i = 1:3
        traceDataInt(:,i) = interp1(tt1, traceData(:,i), tti, 'spline'); 
    end;
    plot3(traceDataInt(:,1), traceDataInt(:,2), ...
          traceDataInt(:,3), '-', 'LineWidth', 1.0, 'Color', lightGray);
else
    plot3(traceData(:,1), traceData(:,2), ...
          traceData(:,3), '-', 'LineWidth', 1.0, 'Color', lightGray);
end;

%% Mark the time samples
for itime = 1:noTime
    color = map(floor(1 + dcolor*timeSample*(itime - 1)), :);
    plot3(traceData(itime,1), traceData(itime,2), traceData(itime,3), 'o', ...
      'LineWidth', 0.5, 'MarkerSize', markerSizeSample, ...
      'MarkerEdgeColor', 'w', 'MarkerFaceColor', color);
end;

%% Mark the beginning and the end of the hodogram
plot3(traceData(1,1), traceData(1,2), traceData(1,3), symbolBeg, ...
      'LineWidth', 0.5, 'MarkerSize', markerSizeBeg, ...
      'MarkerEdgeColor', 'w', 'MarkerFaceColor', map(1,:));
plot3(traceData(end,1), traceData(end,2), traceData(end,3), symbolEnd, ...
      'LineWidth', 0.5, 'MarkerSize', markerSizeEnd, ...
      'MarkerEdgeColor', 'w', 'MarkerFaceColor', map(end,:));


%% Add colorbar to the plot
ch = colorbar('location', 'EastOutside');  
xlabel(ch, ['Time (', timeUnits, ')'], 'Fontsize', 14, 'Fontname', 'Times');
set(gca, 'CLim', timeSample*[timeBeg timeEnd], 'Fontsize', 14, 'Fontname', 'Times')

axis(maxTrace*[-1, 1, -1, 1, -1, 1]);     
axis('vis3d');    view(3);    drawnow;

end    % of the function
