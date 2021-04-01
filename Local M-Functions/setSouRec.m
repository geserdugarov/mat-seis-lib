%% *setSouRec*
% Input locations of sources and receivers

%% 
% *NB:*
%
% * This script needs to be tailored to each specific data set

%%
% *Input:* none

%%
% *Output:*

% xSou      - [3, noSou] array of the source coordinates
% noSou     - [scalar] number of sources
% xRec      - [3, noRec] array of the receiver coordinates
% noRec     - [scalar] number of receivers
% xWell     - [2, 1] array of lateral coordinates of a vertical observation well
% noPerf    - [scalar] number of perforation shots
% noEvnt    - [scalar] number of microseismic events
%             noSou = noPerf + noEvnt
% indexPerf - [1, noSou] array of the indexes of perforation shots, which are found as 
%             find(indexPerf == 1); other elements of indexPerf are 0
% indexEvnt - [1, noSou] array of the indexes of microseismic events, which are found as 
%             find(indexEvnt == 1); other elements of indexEvnt are 0

%%
% *Author:* Vladimir Grechka 2012

%%
function [xSou, noSou, xRec, noRec, xWell, noPerf, noEvnt, indexPerf, indexEvnt] = setSouRec
%% Settings  
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

xWell = [];

global sourecflag;
%% Input the source coordinates
if sourecflag==1
    xSou = [0 0 0]';
elseif sourecflag==2
    xSou = [1000 0 0]';
elseif sourecflag==3
	[sx1,sy1,sz1] = meshgrid(-200:100:200, -200:100:200, 2900:100:3200);
    sx2 = reshape(sx1, 1, numel(sx1));
    sy2 = reshape(sy1, 1, numel(sy1));
    sz2 = reshape(sz1, 1, numel(sz1));
    xSou = [sx2; sy2; sz2];
elseif sourecflag==4 || sourecflag==5 || sourecflag==6 || sourecflag==7
    xSou = [0; 0; 0];
else
    sx = (-1000 : 200 : 1000);  sz = (2800 : 100 : 3200); 
    [sx1, sz1] = meshgrid(sx, sz);
    sx2 = reshape(sx1, 1, numel(sx1));
    sz2 = reshape(sz1, 1, numel(sz1));   sy2 = zeros(size(sx2));  
    xSou   = [sx2; sy2; sz2];
end
noSou  = size(xSou,2);

%% Assign the indexes of perforation shots and events

% One's in indexPerf and indexEvnt correspond to the perforation shots and events;
% the indexes should be complementary
%perfNum = [2];   noPerf = length(perfNum);
perfNum = [];   noPerf = length(perfNum);
indexPerf = zeros(1,noSou);  
indexPerf(1, perfNum) = 1;

noEvnt = noSou - noPerf;   
indexEvnt = ones(1,noSou);
indexEvnt(indexPerf == 1) = 0;

if sum(indexPerf + indexEvnt) ~= noSou
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    disp(indexPerf);
    disp(indexEvnt);
    fprintf('>>> One''s placed in ''indexPerf'' or ''indexEvnt'' are incorrect \n \n');
      error('>>> STOP');
end

%% Input the receiver coordinates
% rx = (-4000 : 1000 : 4000);  ry = rx;  %0;  %rx;  
% [rx1, ry1] = meshgrid(rx, ry);
% rx2 = reshape(rx1, 1, numel(rx1));
% ry2 = reshape(ry1, 1, numel(ry1));   rz2 = zeros(size(rx2));  
% xRec   = [rx2; ry2; rz2];
if sourecflag==1
%     temp = [200:200:2000]; % x-coordinates
%     temp = [400 600];
    temp = 2000;
    xRec = [temp; zeros(1, numel(temp)); zeros(1, numel(temp))];
    clear temp;
elseif sourecflag==2
    temp = 1100:100:1500; % x-coordinates
    xRec = [zeros(1, numel(temp)); zeros(1, numel(temp)); temp];
    clear temp;
elseif sourecflag==4
    temp = 50:50:2000;
    azim = deg2rad(60);
%     temp = 2000;
%     azim = deg2rad(60:45:150);
    xRec = zeros(3,numel(temp)*numel(azim));
    for i = 1:numel(azim)
        xRec(1,(i-1)*numel(temp)+1 : i*numel(temp)) = temp*cos(azim(i));
        xRec(2,(i-1)*numel(temp)+1 : i*numel(temp)) = temp*sin(azim(i));
    end
    clear temp azim i;
    
elseif sourecflag==5 || sourecflag==6
    
    temp = 2000;
    azim = deg2rad(0);
    
    xRec = zeros(3,numel(temp)*numel(azim));
    for i = 1:numel(azim)
        xRec(1,(i-1)*numel(temp)+1 : i*numel(temp)) = temp*cos(azim(i));
        xRec(2,(i-1)*numel(temp)+1 : i*numel(temp)) = temp*sin(azim(i));
    end
    clear temp azim i;
    
elseif sourecflag==7
    temp = 100:100:3500;
    azim = deg2rad(0-90);
%     temp = 1500;
%     azim = deg2rad((0:5:360) - 90);
%     temp = 2500;
%     azim = deg2rad(30-90);
    xRec = zeros(3,numel(temp)*numel(azim));
    for i = 1:numel(azim)
        xRec(1,(i-1)*numel(temp)+1 : i*numel(temp)) = temp*cos(azim(i));
        xRec(2,(i-1)*numel(temp)+1 : i*numel(temp)) = temp*sin(azim(i));
    end
    clear temp azim i;
else
%     temp = 50:50:2000;
%     azim = deg2rad(60:15:150);
%     xRec = zeros(3,numel(temp)*numel(azim));
%     for i = 1:numel(azim)
%         xRec(1,(i-1)*numel(temp)+1 : i*numel(temp)) = temp*cos(azim(i));
%         xRec(2,(i-1)*numel(temp)+1 : i*numel(temp)) = temp*sin(azim(i));
%     end
%     clear temp azim i;
%     noRec = size(xRec,2);
end
noRec = size(xRec,2);

end    % of the function
