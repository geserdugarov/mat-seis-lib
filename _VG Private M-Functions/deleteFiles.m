%% *deleteFiles*
% Delete files with given prefix and suffix from a folder and all its subfolders

%% 
% *Input:*

% rootFolder - [string] full path to the folder from which files whose names begin with 
%              'filePrefix' and end with 'fileSuffix' are to be deleted 
% filePrefix - [string] prefix of files 
% fileSuffix - [string] suffix of files 

%%
% *Output:* none

%%
% *Author:* Vladimir Grechka 2012 - 2014

%%
function deleteFiles(rootFolder, filePrefix, fileSuffix)
%% Settings  
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

recycle on;    % 'recycle on' makes Matlab 'delete' function place the deleted files in the 
               % 'Recycle Bin' rather than delete them permanently

%% Find files to be deleted
[fileList, fileCount] = findFiles(rootFolder, filePrefix, fileSuffix);
display(fileList);

if fileCount ~= 0
    beep on;  beep;  beep off;
    fprintf('>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    tmp = input('>>> Please type ''yes'' to confirm the deletion of the above listed file(s) -> ');  

    if strcmp(tmp, 'yes') == 1
        display('>>> Working. Please wait...');
        errorFlag = 1;  
        while errorFlag == 1;  
            try 
                for ifile = 1:fileCount
                    delete(fileList{ifile});
                end;
                fprintf('>>> %g file(s) from folder ''%s'' have been deleted \n', ...
                    fileCount, rootFolder);   
                errorFlag = 0; 
            catch err; 
                errorFlag = 1;  
                fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
                display('>>> Error message:')
                fprintf('>>> %s \n', err.message);
            end; 
        end;
    else
        fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
        display('>>> No file(s) will be deleted');  beep on;  beep;  beep off;  
    end;
    
else
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    display('>>> Found no file(s) to delete');  beep on;  beep;  beep off;  
end;

end    % of the function

%%
