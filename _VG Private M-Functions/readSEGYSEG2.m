%% *readSEGYSEG2*
% Read seismic traces in SEGY or SEG2 format and save them as .mat file 

%%
% *Input:*

% fullFileName        - [string] full name of file containing 3C seismic data in .seg2 
%                       or .segy format
% flagComp            - [string] equal to '1C' or '3C' and denoting the number of components in 
%                       the input file
% flagSave            - flag equal to 1 to save traces as .mat file and ~= 1 otherwise    
% outputFolder        - [string] full name of the folder in which traces in .mat are to be saved
%                       (*) if (flagSave == 1 && isempty(outputFolder) == 1), '.mat' file will be 
%                           written in directory 'Outputs' under the directory from which the main
%                           script is invoked
% outputTimeUnits     - [string] units of the output time sample
% flagNaN             - [scalar] flag equal to 1 to replace NaNs in traces with zeros and ~= 1 
%                       to leave them as NaNs    

%%
% *Output:*

% traceData           - [header.noSample, header.noTrace] trace array
%
% header              - structure containing the following header information:        
%   header.fileName   - [string] name of seismic file
%   header.filePath   - [string] full name of the folder containing the file 
%   header.noTrace    - [scalar] the number of traces
%   header.noRec      - [scalar] the number of receivers (= header.noTrace/3)
%   header.noSample   - [scalar] the number of time samples   
%   header.timeSample - [string] time sample in outputTimeUnits 
%                       (*) if (isempty(outputTimeUnits) == 1), 'header.timeSample' is in ms
%   header.timeUnits  - [string] units of the time sample  
%
%   header.acquisition.date           - [scalar] date of data acquisition whose numeric value
%                                       is eqaul to day/month/year
%   header.acquisition.time           - [string] time of data acquisition in format 'HH:MM:SS'
%   header.acquisition.secondFraction - [scalar] fraction of a second of the first sample in input
%                                       file

%%
% *Author:* Vladimir Grechka 2012 - 2014

%%
function [traceData, header] = ...
    readSEGYSEG2(fullFileName, flagComp, flagSave, outputFolder, outputTimeUnits, flagNaN)
%% Settings and defaults
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
narginTrue = nargin(thisFileName);   
narginchk(narginTrue, narginTrue);

header.noSample = [];   header.timeSample = [];
header.noTrace  = [];   header.noRec      = [];   header.timeUnits  = 'ms';
header.acquisition.date = [];    header.acquisition.time = [];    
header.acquisition.secondFraction = [];
traceData = [];   msg1 = [];   msg2 = [];
                   
if isempty(flagComp) == 1;
    flagComp = '3C';
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Parameter ''flagComp'' is defaulted to ''%s'' \n', flagComp);
end;    

global letterDrive

%% Break 'fullFileName' into the folder, file name, and extension
[seismicFolderName, seismicFileName, seismicFileExtention] = fileparts(fullFileName);
header.fileName =  [seismicFileName, seismicFileExtention];
header.filePath =   seismicFolderName;

seismicFormat = [0, 0];
seismicFormat(1) = ...
    strcmp(seismicFileExtention, '.sg2') + strcmp(seismicFileExtention, '.SG2') + ...
    strcmp(seismicFileExtention, '.dat') + strcmp(seismicFileExtention, '.DAT') + ...
    strcmp(seismicFileExtention, '.seg2');
seismicFormat(2) = ...
    strcmp(seismicFileExtention, '.sgy') + strcmp(seismicFileExtention, '.segy') + ...
    strcmp(seismicFileExtention, '.pset');

%% Read .seg2 or .segy file
if seismicFormat(1) > 0
    %% SEG2
    try                                                     % reading with 'Seg2FileReader'
        handleSet   = Seg2FileReader(fullFileName);
        traceHeader = readTraceHeader(handleSet);
        traceData   = readTraceData(handleSet);

%         header.acquisition.date = handleSet.FileHeader.TextData.ACQUISITION_DATE;    
%         header.acquisition.time = handleSet.FileHeader.TextData.ACQUISITION_TIME;    
%         header.acquisition.secondFraction = handleSet.FileHeader.TextData.ACQUISITION_SECOND_FRACTION;
        
        if ischar(traceHeader(1,1).TextData.SAMPLE_INTERVAL) == 1
            i1 = find(uint8(traceHeader(1,1).TextData.SAMPLE_INTERVAL) == 1);
            header.timeSample = 1.e+3*str2double( ...
                                      traceHeader(1,1).TextData.SAMPLE_INTERVAL(1:i1(1)-1));     
        else
            header.timeSample = 1.e+3*traceHeader(1,1).TextData.SAMPLE_INTERVAL;  
        end;

    catch ME  
        msg1  = ME.message;
        line1 = ME.stack.line;
    end;

    if isempty(msg1) == 0
       fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
       fprintf('>>> Failure of ''Seg2FileReader'' to read file \n');
       fprintf('    ''%s'' \n', fullFileName);
       fprintf('>>> ERROR on line %g: %s \n', [line1, msg1]);

        clear ME;
        try
            [traceData, headerSEG2] = seg2load(fullFileName);   % reading with 'seg2load'
            if isempty(header.timeSample) == 1
                header.timeSample = 1.e+3*headerSEG2.tr.sampling(1);
            end;
        catch ME
            msg2  = ME.message;
            line2 = ME.stack.line;
        end;
    end;
          
    if isempty(msg2) == 0
        fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
        fprintf('>>> Failure of ''seg2load'' to read file \n');
        fprintf('    ''%s'' \n', fullFileName);
        fprintf('>>> ERROR on line %g: %s \n', [line2, msg2]);
    end;
    
elseif seismicFormat(2) > 0
    %% SEGY
    try                                                     % reading with 'SegYFileReader'
        handleSet   = SegYFileReader(fullFileName);
        traceHeader = readTraceHeader(handleSet);
        traceData   = readTraceData(handleSet);
        header.timeSample = 1.e-3*traceHeader(1,1).SampleInterval_ms;
    catch ME  
        msg1  = ME.message;
        line1 = ME.stack.line;
    end;
    
    if isempty(msg1) == 0
        fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
        fprintf('>>> Failure of ''SegYFileReader'' to read file \n');
        fprintf('    ''%s'' \n', fullFileName);
        fprintf('>>> ERROR on line %g: %s \n', [line1, msg1]);
        
        clear ME;
        try                                                 % reading with 'ReadSegy'
            localFolderName = fullfile(letterDrive, 'Playgrounds\SegyMAT');
            addpath(localFolderName);
            [traceData, ~, SegyHeader] = ReadSegy(fullFileName);          
            header.timeSample = 1.e-3*SegyHeader.dt;   % time sample in ms
            rmpath(localFolderName);          
        catch ME
            msg2  = ME.message;
            line2 = ME.stack.line;
        end;
    end;
    
    if isempty(msg2) == 0
        fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
        fprintf('>>> Failure of ''ReadSegy'' to read file \n');
        fprintf('    ''%s'' \n', fullFileName);
        fprintf('>>> ERROR on line %g: %s \n', [line2, msg2]);
    end;

else
    %% Catch the error
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Incorrect extension of seismic data files \n');
    fprintf(['>>> Legitimate extensions are: ''.sg2'', ''.SG2'', ''.seg2'', ''.dat'', ', ...
             '''.segy'', ''.sgy'', ''.pset'' \n \n']);
      error('>>> STOP');
end;

if isempty(msg1) == 0 && isempty(msg2) == 0
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Cannot read file \n');
    fprintf('    ''%s'' \n', fullFileName);
%    fprintf('>>> PAUSE -- Continue? \n');   pause;
end;
   
%% Fill the output header
header.noSample = size(traceData, 1);
header.noTrace  = size(traceData, 2);
if strcmp(flagComp, '3C') == 1
    header.noRec = header.noTrace/3;    % 3C data
else
    header.noRec = header.noTrace;      % 1C data
end;

%% Change the time units
if isempty(outputTimeUnits) == 0
    header.timeSample = u2u(header.timeSample, 'ms', outputTimeUnits, []);
    header.timeUnits  = outputTimeUnits;
end;

%% Replace NaN trace values with zeros
if flagNaN == 1;   
    traceData(isnan(traceData)) = 0;   
end;

%% Output .mat structure   
if flagSave == 1
    if isempty(outputFolder) == 1
        outputFolder = fullfile(pwd, 'Outputs');
    end;
    if exist(outputFolder, 'dir') == 0
        fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
        fprintf('>>> Creating folder ''%s'' \n', outputFolder);
        fprintf('>>>    to save file ''%s'' \n \n', [seismicFileName, '.mat']);
        mkdir(outputFolder);
    end;
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Saving file ''%s'' \n', [seismicFileName, '.mat']);
    fprintf('    in folder ''%s'' \n', outputFolder);
    save(fullfile(outputFolder, [seismicFileName, '.mat']), '-mat', 'traceData', 'header');
end;

end    % of the function

