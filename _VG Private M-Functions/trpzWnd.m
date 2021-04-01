%% *trpzWnd*
% Generate a trapezoidal window 

%%
% *Input:*

% sbeg - [integer] number of the first time sample of created window
% send - [integer] number of the last time sample of created window
% st   - [1, 4][integer] array of samples defining the abscissae of a trapezoidal
%        (*) It is assumed that st(i) are sorted in the ascending order 
% yt   - [1, 2] array of ordinates of the parallel sides of a trapezoidal
%        (*) It is assumed that yt(1) < yt(2) 

%%
% *Output:*

% trpz - [1, send-sbeg] array of the trapezoidal ordinates

%%
% *Author:* Vladimir Grechka 2012

%%
function [trpz] = trpzWnd(sbeg, send, st, yt)
%% Seetings and checks
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

% Check whether st and yt are properly sorted
ss = sort(st);
if isequal(st, ss) == 0
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Array st = [%g, %g, %g, %g] is not sorted in the ascending order \n \n', st);
      error('>>> STOP');
end;

if yt(2) < yt(1)
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> The condition yt(2) > yt(1) for yt = [%g, %g] is violated \n \n', yt);
      error('>>> STOP');
end;

% Various checks
if length(st) ~= 4
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Array ''st'' has the length %g instead of 4 \n \n', length(st));
      error('>>> STOP');
end;

if length(yt) ~= 2
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Array ''yt'' has the length %g instead of 2 \n \n', length(yt));
      error('>>> STOP');
end;

if st(1) < sbeg
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Trapezoidal is outside the sample range: st(1) = %g, sbeg = %g \n \n', ...
            [st(1), sbeg]);
      error('>>> STOP');
end;

if st(4) > send
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Trapezoidal is outside the sample range: st(4) = %g, send = %g \n', ...
            [st(4), send]);
      error('>>> STOP');
end;

%% Generate a trapezoidal window
s = sbeg:1:send;
trpz = yt(1)*ones(1, length(s));                    % horizontal sides
trpz(st(2) <= s  &  s <= st(3)) = yt(2);

i1 = find(st(1) < s  &  s < st(2));                 % dipping left side
trpz(i1) = yt(1) + (yt(2) - yt(1))/(st(2) - st(1))*(s(i1) - st(1));

i2 = find(st(3) < s  &  s < st(4));                 % dipping right side
trpz(i2) = yt(1) + (yt(2) - yt(1))/(st(3) - st(4))*(s(i2) - st(4));

end    % of the function