%% *cell2dateNumber*
% Covert a cell array containing date information to a column of date numbers

%%
% *Input:* 

%   cellArray  - [cell] array containing date information
%   pointer    - [1, 3] array positions in string 'cellArray{i}', at which the date string
%                (in 'yyyymmdd' format), the time string (in 'hhmnss' format), and the  
%                millisecond (in 'msc' format) start  

%%
% *Output:*

%   dateVector - [size(cellArray, 1), 6] array, whose rows contain numeric elements
%                [yyyy, mm, dd, hh, mn, ss + 0.001*msc]
%   dateNumber - [size(cellArray, 1), 1] vector of the date numbers given by Matlab function
%                'datenum'

%%
% *Author:* Vladimir Grechka 2013

%%
function [dateVector, dateNumber] = cell2dateNumber(cellArray, pointer)
%% Settings 
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);
    
dateVector = NaN(size(cellArray, 1), 6);
dateNumber = NaN(size(cellArray, 1), 1);

%% Checks
if isempty(cellArray) == 1
    dateNumber = [];
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Empty input ''cellArray'' -> empty output ''dateNumber'' --> PAUSE \n');  pause;
end;

%% Construct the date-number array
for i = 1:size(cellArray, 1)
    %% Date
    if isnan(pointer(1)) == 1
        fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
        fprintf('>>> WARNING: Undefined date ''yyyymmdd'' --> set to zero \n');  pause;
        yyyy = 0;  mm = 0;  dd = 0;
    else
        if pointer(1) + 7 > size(cellArray{i}, 2)
            fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
            fprintf(['>>> size(cellArray{%g}, 2) = %g is to small to accommodate ''yyyymmdd'' ', ...
                     'starting at position number %g \n'], i, size(cellArray{i}, 2), pointer(1));  
              error('>>> STOP');
        else
            yyyy = str2double(cellArray{i}(pointer(1)     : pointer(1) + 3));
            if isnan(yyyy) == 1
                fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
                fprintf('>>> Varible ''yyyy'' (= %s) extracted from ''%s'' \n', ... 
                    cellArray{i}(pointer(1) : pointer(1) + 3), cellArray{i} ); 
                fprintf('    is not a number \n');
                 error('>>> STOP');
            end;
            mm = str2double(cellArray{i}(pointer(1) + 4 : pointer(1) + 5));
            if isnan(mm) == 1
                fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
                fprintf('>>> Varible ''mm'' (= %s) extracted from ''%s'' \n', ... 
                    cellArray{i}(pointer(1) + 4 : pointer(1) + 5), cellArray{i} ); 
                fprintf('    is not a number \n');
                 error('>>> STOP');
            end;
            dd = str2double(cellArray{i}(pointer(1) + 6 : pointer(1) + 7));
            if isnan(dd) == 1
                fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
                fprintf('>>> Varible ''dd'' (= %s) extracted from ''%s'' \n', ... 
                    cellArray{i}(pointer(1) + 6 : pointer(1) + 7), cellArray{i} ); 
                fprintf('    is not a number \n');
                 error('>>> STOP');
            end;
        end;
    end;
    
    %% Time
    if isnan(pointer(2)) == 1
        fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
        fprintf('>>> WARNING: Undefined time ''hhmnss'' is set to zero \n');   
        hh = 0;  mn = 0;  ss = 0;
    else
        if pointer(2) + 5 > size(cellArray{i}, 2)
            fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
            fprintf(['>>> size(cellArray{%g}, 2) = %g is to small to accommodate ''hhmnss'' ', ...
                     'starting at position number %g \n'], i, size(cellArray{i}, 2), pointer(2));  
              error('>>> STOP');
        else
            hh = str2double(cellArray{i}(pointer(2)     : pointer(2) + 1));
            if isnan(hh) == 1
                fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
                fprintf('>>> Varible ''hh'' (= %s) extracted from ''%s'' \n', ... 
                    cellArray{i}(pointer(2) : pointer(2) + 1), cellArray{i} ); 
                fprintf('    is not a number \n');
                 error('>>> STOP');
            end;
            mn = str2double(cellArray{i}(pointer(2) + 2 : pointer(2) + 3));
            if isnan(mn) == 1
                fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
                fprintf('>>> Varible ''mn'' (= %s) extracted from ''%s'' \n', ... 
                    cellArray{i}(pointer(2) + 2 : pointer(2) + 3), cellArray{i} ); 
                fprintf('    is not a number \n');
                 error('>>> STOP');
            end;
            ss = str2double(cellArray{i}(pointer(2) + 4 : pointer(2) + 5));
            if isnan(ss) == 1
                fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
                fprintf('>>> Varible ''ss'' (= %s) extracted from ''%s'' \n', ... 
                    cellArray{i}(pointer(2) + 4 : pointer(2) + 5), cellArray{i} ); 
                fprintf('    is not a number \n');
                 error('>>> STOP');
            end;
        end;
    end;

    %% Millisecond
    if isnan(pointer(3)) == 1
        msc = 0;   
    else
        if pointer(3) + 2 > size(cellArray{i}, 2)
            msc = 0;
        else
            msc = str2double(cellArray{i}(pointer(3) : pointer(3) + 2));
            if isnan(msc) == 1
                fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
                fprintf('>>> WARNING: Varible ''msc'' (= %s) extracted from ''%s'' \n', ... 
                    cellArray{i}(pointer(3) : pointer(3) + 2), cellArray{i} ); 
                fprintf('             is not a number \n');
                 error('>>> STOP');
            end;
        end;
    end;

    %% The dateNumber element
    dateVector(i, :) = [yyyy, mm, dd, hh, mn, ss + 0.001*msc];
    dateNumber(i, 1) = datenum(dateVector(i,:));
        
end;

end  % of the function

%%

