%% *positiveEnvelope*
% Given a cloud y(x), find the x-y pairs that envelope the cloud from above

%%
% *Input:*

% x        - array of abscissas of points that form cloud y(x)
% y        - array of ordinates of points that form cloud y(x)
% relWidth - [scalar] determining the relative width of the envelope (0 <= relWidth <= 1)

%%
% *Output:*

% xenv     - array of abscissas of the envelope
% yenv     - array of ordinates of the envelope

%%
% *Author:* Vladimir Grechka 2013

%% 
% *Acknowledgements:* Modified Nicolas Hummel's algorithm (see 'triggering_front_fitting.m')

%%
function [xenv, yenv] = positiveEnvelope(x, y, relWidth)
%% Settings and checks
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

if relWidth < 0  ||  relWidth > 1
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Improper value of parameter relWidth = %g \n', relWidth);
    fprintf('    Its legitimate values should satisfy the inequality 0 <= relWidth <= 1 \n');                
      error('>>> STOP');
end;  

%% Sort the data
[xx, ix] = sort(x, 'ascend'); 
yy = y(ix);

%% Find the envelope
count = 1;  
xenv(1,1) = xx(1);  
yenv(1,1) = yy(1);  

for i = 2:length(xx)
    if yy(i) > (1 - relWidth)*max(yenv)
        count = count + 1;  
        xenv(1,count) = xx(i);  
        yenv(1,count) = yy(i);
    end;
end;

end  % of the function