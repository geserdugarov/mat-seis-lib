%% *u2u*
% Unit-to-unit converter

%% 
% *Input:*

% quantityInp - [scalar, vector, or matrix] quantity whose units are to be converted   
% unitsInp    - [char] input units  
% unitsOut    - [char] output units 
% warningFlag - [scalar] equal to 0 or 1 that turns off and on the warning message pertaining
%               to the output NaN values

%%
% *Output:* 

% quantityOut - [scalar, vector, or matrix] quantityInp in output units  

%%
% *Supported units:*

% length    - m, km, ft, kft 
% time      - s, ms, us, min, h
% frequency - Hz, kHz
% angle     - deg, rad
% mass      - g, kg, lb

%%
% *Author:* Vladimir Grechka 2012

%%
function [quantityOut] = u2u(quantityInp, unitsInp, unitsOut, warningFlag)
%% Settings and defaults
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
% narginTrue = nargin(thisFileName);    % There two lines are commented out because typing the 
% narginchk(narginTrue, narginTrue);    % 'warningFlag' in each 'u2u' call gets annoying

quantityOut = quantityInp;

% Cell aray of the supported units  
suppUnits = {  'm',  'km', 'ft', 'kft', ...         % length
               's',  'ms', 'us', 'min', 'h', ...    % time 
              'Hz', 'kHz', ...                      % frequency
             'deg', 'rad', ...                      % angle
               'g',  'kg', 'lb'};                   % mass

if nargin == 3;    warningFlag = 0;    end;
    
%% Check whether the input and output units exist
if isempty(unitsInp) == 1 
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Variable ''unitsInp'' is undefined --> the input units remain unchanged \n \n');    
    return;
end;

if isempty(unitsOut) == 1 
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Variable ''unitsOut'' is undefined --> the input units remain unchanged \n \n');
    return;
end;

%% Check whether the input and output units are supported 
if ischar(unitsInp) == 0 || ismember({unitsInp}, suppUnits) == 0
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Variable unitsInp = ''%s'' \n', unitsInp);
    fprintf('>>> is either not a character or currently unsupported \n \n');
      error('>>> STOP');
end; 

if ischar(unitsOut) == 0 || ismember({unitsOut}, suppUnits) == 0
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Variable unitsOut = ''%s'' \n', unitsOut);
    fprintf('>>> is either not a character or currently unsupported \n \n');
      error('>>> STOP');
end; 

%% Create the unit-conversion chart
f = diag(ones(size(suppUnits,2),1));    % place 1 on the diagonal of matrix f
f(f == 0) = NaN;                        % fill the rest of matrix f with NaNs

% Fill the upper-right corner of 'f' with the conversion factors of units in the order 
% given be 'suppUnits':
% +----+----+----+-----+---+----+----+-----+---+----+-----+-----+-----+----+----+----+
% |  1 |  2 |  3 |  4  | 5 |  6 |  7 |  8  | 9 | 10 |  11 |  12 |  13 | 14 | 15 | 16 |
% |  m | km | ft | kft | s | ms | us | min | h | Hz | kHz | deg | rad |  g | kg | lb |
% +----+----+----+-----+---+----+----+-----+---+----+-----+-----+-----+----+----+----+

% Length
m2ft = 3.28084;    
f(1,2) = 0.001;             % m2km
f(1,3) = m2ft;              % m2ft
f(1,4)  = 0.001*m2ft;       % m2kft
f(2,3) = 1000*m2ft;         % km2ft
f(2,4)  = m2ft;             % km2kft
f(3,4)  = 0.001;            % ft2kft    

% Time
f(5,6) = 1000;              % s2ms
f(5,7) = 1.e+6;             % s2us
f(5,8) = 1/60;              % s2min
f(5,9) = 1/3600;            % s2h
f(6,7) = 1000;              % ms2us
f(6,8) = 1.e-3/60;          % ms2min
f(6,9) = 1.e-3/3600;        % ms2h
f(7,8) = 1.e-6/60;          % us2min
f(7,9) = 1.e-6/3600;        % us2h
f(8,9) = 1/60;              % min2h

% Frequency
f(10,11) = 0.001;             % Hz2kHz

% Angle
f(12,13) = pi/180;          % deg2rad  

% Mass
kg2lb = 2.204623;           
f(14,15) = 0.001;           % g2kg
f(14,16) = 0.001*kg2lb;     % g2lb  
f(15,16) = kg2lb;           % kg2lb

%% Fill the symmetric part with reciprocal factors
for i=1:size(f,2)
    for j=1:i-1
        f(i,j) = 1/f(j,i);
    end;
end;

%% Convert the input quantities
[~, jInp] = ismember({unitsInp}, suppUnits);
[~, jOut] = ismember({unitsOut}, suppUnits);
quantityOut = f(jInp, jOut)*quantityInp;

%% Check whether any NaNs have been assigned
iNaN = find(isnan(reshape(quantityOut, numel(quantityOut), 1)));
if sum(iNaN) ~= 0  && warningFlag == 1
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> WARNING: variable ''quantityOut'' contains NaN values \n');
    display(quantityOut);
end;
%if warningFlag == 1;    display(quantityOut);    end;

end    % of the function
