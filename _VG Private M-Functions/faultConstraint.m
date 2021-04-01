%% *faultConstraint*
% Create the linear constraint inequalty for a fault reflection point 

%%
% *Input:*

% xy       - [2*(noSeg-1), 1] array containing pairs of the x- and y-coordinates of 
%            intersections of the ray trajectory with model interfaces and faults
% zInt     - [3, :] z-coefficients of interfaces produced by 'planeEquation'
% fltLayer - [scalar] layer number in which the fault reflection takes place 
% rayCode  - [3, noSeg] array resulting from 'produceRayCode'
%            noSeg - the number of ray segments

%%
% *Output:*

% Acon     - [4*(noSeg-1), 2*(noSeg-1)] rectangular matrix of linear constraints 
%            dot(Acon, xy) < Bcon as required in Matlab function 'fmincon' 
%            The top and bottom blocks of Acon take care of the top and bottom interfaces 
%            bounding a fault segment 
% Bcon     - [4*(noSeg-1), 1] vector of linear constraints dot(Acon, xy) < Bcon as required 
%            in Matlab function 'fmincon' 

%%
% *Author:* Vladimir Grechka 2012 

%%
function [Acon, Bcon] = faultConstraint(xy, zInt, zFlt, fltLayer, rayCode)
%% Settings  
[~, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

%%
noInt = size(zInt, 2);              % the number of model interfaces   
xyLength = size(xy, 1);
Acon = zeros(2*xyLength, xyLength);   
Bcon =  ones(2*xyLength, 1);

iseg = find(rayCode(4,:) ~= 0);     % ray segment at which reflection from a fault occurs
numbrFlt = rayCode(4,iseg);         % the fault number

%%
if fltLayer ~= 1 
    % Inequality for the top interface is invoked; top blocks of Acon and Bcon are in action
    zTopInt = dot(zInt(1:2,fltLayer-1), xy(2*iseg - 1 : 2*iseg, 1)) + zInt(3,fltLayer-1);
    Acon(iseg, 2*iseg - 1 : 2*iseg) = -zFlt(1:2,numbrFlt);
    Bcon(iseg, 1) = zFlt(3,numbrFlt) - zTopInt;
end;

%%
if fltLayer ~= noInt+1      
    % Inequality for the bottom interface is invoked; bottom blocks of Acon and Bcon are in action
    zBotInt = dot(zInt(1:2,fltLayer), xy([2*iseg - 1 : 2*iseg], 1)) + zInt(3,fltLayer);
    Acon(xyLength + iseg, [2*iseg - 1 : 2*iseg]) = zFlt(1:2,numbrFlt);
    Bcon(xyLength + iseg, 1) = zBotInt - zFlt(3,numbrFlt);
end;

end    % of the function