%% *findFiles*
% List files with given prefix and suffix in folder and all its subfolders

%% 
% *Input:*

% rootFolder - [string] full path to the folder in which file names beginning with 'filePrefix' 
%              and ending with 'fileSuffix' are to be searched 
% filePrefix - [string] prefix of files 
% fileSuffix - [string] suffix of files 

%%
% *Output:*

% fileList   - [cell] list of full names of files in 'rootFolder' and all its subfolders 
%              beginning with 'filePrefix' and ending with 'fileSuffix'
% fileCount  - [scalar] number of files in 'fileList'

%%
% *Author:* Vladimir Grechka 2012 2013

function [fileList, fileCount] = findFiles(rootFolder, filePrefix, fileSuffix)
%% Settings and defaults 
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

fileList  = {};   fileCount = 0;

if isempty(fileSuffix) == 1;
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Empty variable ''fileSuffix'' \n');   
end;

if isempty(filePrefix) == 1;   
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Empty variable ''filePrefix'' \n');  
end;

fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
fprintf('>>> Searching for files with prefix ''%s'' and suffix ''%s'' in folder \n', ...
             filePrefix, fileSuffix);
fprintf('    ''%s'' \n', rootFolder);

%% Form a string that includes all folders and subfolders below 'rootFolder'
stringOfFolders = genpath(rootFolder);

% Break 'stringOfFolders' into strings containing individual subfolders
indexSemiColon = find(stringOfFolders == ';');

if isempty(indexSemiColon) == 1
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Folder \n');
    fprintf('    ''%s'' \n', rootFolder);
    fprintf('    is not found \n \n');
      error('>>> STOP');
end;

subFolder(1,1) = {stringOfFolders(1 : indexSemiColon(1)-1)};
if length(indexSemiColon) > 1
    for i = 2:length(indexSemiColon)
        subFolder(i,1) = {stringOfFolders(indexSemiColon(i-1)+1 : indexSemiColon(i)-1)};
    end;
end;

%% Search subfolders for file names beginning with 'filePrefix' ending with 'fileSuffix'
for i = 1:size(subFolder,1) 
    clear listOfFiles fileName
    listOfFiles = dir(subFolder{i,1});
    fileName = {listOfFiles.name}';
    
    % Search the i-th subfolder for files
    for j=1:size(fileName,1)

        if isempty(filePrefix) == 0  &&  isempty(fileSuffix) == 0
            % Searching files by both suffix and prefix
            if length(fileName{j,1}) >= length(filePrefix) + length(fileSuffix)
                if strcmp(fileName{j,1}(1 : length(filePrefix)), filePrefix) == 1  && ...
                   strcmp(fileName{j,1}(end+1-length(fileSuffix) : end), fileSuffix) == 1
                    fileCount = fileCount + 1;
                    fileList(fileCount,1) = {fullfile(subFolder{i,1}, fileName{j,1})};
                end;
            end;
        end;
        
        if isempty(filePrefix) == 0  &&  isempty(fileSuffix) == 1
            % Searching files by prefix only
            if length(fileName{j,1}) >= length(filePrefix)
                if strcmp(fileName{j,1}(1 : length(filePrefix)), filePrefix) == 1
                    fileCount = fileCount + 1;
                    fileList(fileCount,1) = {fullfile(subFolder{i,1}, fileName{j,1})};
                end;
            end;
        end;
        
        if isempty(filePrefix) == 1  &&  isempty(fileSuffix) == 0
            % Searching files by suffix only
            if length(fileName{j,1}) >= length(fileSuffix)
                if strcmp(fileName{j,1}(end+1-length(fileSuffix) : end), fileSuffix) == 1
                    fileCount = fileCount + 1;
                    fileList(fileCount,1) = {fullfile(subFolder{i,1}, fileName{j,1})};
                end;
            end;
        end;
        
    end;    % of the file loop
end;    % of the subfolder loop

fprintf('\n>>> %g file(s) with prefix ''%s'' and suffix ''%s'' found \n \n', ...
        fileCount, filePrefix, fileSuffix);
end    % of the function
