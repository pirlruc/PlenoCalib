classdef LytroDesktop
    %LYTRODESKTOP
    %   Lytro Desktop utilities. This assumes that lytro image files are
    %   given in separate folders.
    
    properties
        folder = pwd            % Lytro images folder path
    end
    
    properties (Dependent)
        folders                 % Subfolders in folder path
        numberFolders           % Number of subfolders in folder path
        files                   % Lightfield files in folder path (.lfp, .raw)
        numberFiles             % Number of lightfield files in folder and subfolders
        metadata                % Lightfield files metadata in folder path
    end
    
    methods
        function self = LytroDesktop(varargin)
            %
            % Create lytro desktop instance.
            % 
            % INPUTS:
            %   1. folderPath - lytro images folder path.
            %
            narginchk(0,1);
            
            if ~isempty(varargin)
                if nargin >= 1
                    self.folder = varargin{1};
                end
            end
        end
        
        function self = set.folder(self,newFolderPath)
            self.folder = newFolderPath;
        end
        
        function files = get.files(self)
            % Obtain lightfield image files with extension .lfr, .lfp and 
            % .raw
            files_lfr = dir([ self.folder filesep '**' filesep '*.lfr' ]);
            files_lfp = dir([ self.folder filesep '**' filesep '*.lfp' ]);
            files_raw = dir([ self.folder filesep '**' filesep '*.raw' ]);
            
            % If there are no lightfield files, create empty list
            if size(files_lfr,1) == 0 && size(files_lfp,1) == 0 && size(files_raw,1) == 0
                files = [];
            else
                files = [files_lfr; files_lfp; files_raw];
            end
        end
        
        function number = get.numberFiles(self)
            number = length(self.files);
        end
        
        function folders = get.folders(self)
            % Obtain list of files and folders in directory
            foldersAndFiles = dir([ self.folder filesep '**' filesep '*' ]);
            
            % Obtain only the folders
            isDirectory = [foldersAndFiles.isdir];
            folders     = foldersAndFiles(isDirectory);
            
            % Remove navigation folders
            removeIndexes = cellfun( @strcmp, {folders.name} ...
                                            , repmat({'.'} , 1, length(folders) ));
            folders = folders(~removeIndexes);
            removeIndexes = cellfun( @strcmp, {folders.name} ...
                                            , repmat({'..'}, 1, length(folders) ));
            
            % If there are no ligthfield folder files, create empty list
            if all(removeIndexes)
                folders = [];
            else
                folders = folders(~removeIndexes);
            end
        end
        
        function number = get.numberFolders(self)
            number = length(self.folders);
        end
        
        function metadata = get.metadata(self)
            % Obtain lightfield image files
            metadata = self.files;

            for iFile = 1:self.numberFiles
                % Obtain lightfield filename and .json filename for
                % metadata
                [~,filename]  = fileparts(metadata(iFile).name);
                filepath_json = [metadata(iFile).folder filesep filename '.json'];
                
                % Read metadata from .json file, otherwise read it from
                % lightfield file
                if exist(filepath_json,'file')
                    lytroData.Metadata = LFReadMetadata(filepath_json);
                else
                    % Obtain lightfield path
                    filepath = [metadata(iFile).folder filesep metadata(iFile).name];
                
                    % Obtain filename and extension
                    [~,filename,extension] = fileparts(metadata(iFile).name);
                    
                    % Read metadata associated with raw image file
                    if strcmpi(extension,'.raw') > 0
                        % Obtain metadata file
                        metadataSuffix     = '_metadata.json';
                        metadataFilename   = [metadata(iFile).folder filesep filename metadataSuffix];
                        fileMetadata       = LFReadMetadata(metadataFilename);
                        lytroData.Metadata = fileMetadata;
                    % Read lytro image file
                    else
                        lytroData = LFReadLFP(filepath);
                    end
                end
                
                % Add metadata to file structure
                metadata(iFile).metadata = lytroData.Metadata;
            end
        end
    end
    
    methods
        function deleteUnnecessaryFiles(self)
            %
            % Delete stack view file in lytro image files directory.
            %
            FILENAME_TO_DELETE = ['**' filesep 'stack.view.lfp'];

            % Obtain list of files and folders in current directory
            listFiles = dir([self.folder filesep FILENAME_TO_DELETE]);

            for iFile = 1:length(listFiles)
                % Obtain filepath to delete
                filepath = [listFiles(iFile).folder filesep listFiles(iFile).name];

                % Delete stack view file
                delete(filepath);
            end
        end
        
        function deleteEmptyFolders(self)
            %
            % Delete empty folders in lytro image files directory.
            %

            % Obtain subfolders in current directory
            subfolders = self.folders;
            
            for iFolder = 1:self.numberFolders
                % Obtain folder path
                folderPath = [subfolders(iFolder).folder filesep subfolders(iFolder).name];

                % Remove directory. Matlab only deletes folder if it is
                % empty.
                try
                    rmdir(folderPath);
                catch
                    continue
                end
            end
        end
        
        function sortLightfieldFiles(self)
            %
            % Sort lytro images according to acquisition timestamp and
            % rename folders according to a counter identifier.
            %
            
            % Obtain metadata for image files
            lytroFiles = self.metadata;

            % Obtain acquisition datetime for each lightfield image file
            datetimes  = repmat(datetime,self.numberFiles,1);
            for iFile = 1:self.numberFiles
                % Obtain acquisition datetime for lightfield image
                lytroDatetime = lytroFiles(iFile).metadata.devices.clock.zuluTime;

                % Add datetime. Do not consider miliseconds and time zone.
                datetimes(iFile) = datetime( lytroDatetime(1:end - 5) ...
                                           , 'InputFormat', 'uuuu-MM-dd''T''HH:mm:ss');
            end

            % Sort datetimes of lightfield image files
            [~,orderedIndexes] = sort(datetimes);

            % Rename folders that include lightfield images by a count 
            % identifier. Remember that we have one folder for one
            % ligthfield image.
            for iFile = 1:self.numberFiles
                % Obtain subfolder name of current lightfield image
                [rootFolder,folderName] = fileparts(lytroFiles(orderedIndexes(iFile)).folder);
                newFolderNumber         = sprintf('%05d',iFile);
                
                % If folder name corresponds to the new folder name with
                % counter do not change it
                if strcmp(folderName,newFolderNumber) == true
                    continue
                    
                % Otherwise, rename folder with current counter.
                else
                    folderPath    = lytroFiles(orderedIndexes(iFile)).folder;
                    newFolderPath = [rootFolder filesep newFolderNumber];
                    movefile(folderPath, newFolderPath);
                end
            end
        end
        
        function renameLightfieldFiles(self,filenames)
            %
            % Rename folders of lightfield files.
            %
            % INPUTS:
            %   1. filenames - new folder names for lytro images
            %   subfolders.
            %
            narginchk(2,2);
            
            % Obtain list of lightfield image files
            lytroFiles = self.files;

            % Validate if number of subfolder names are the same of the 
            % number of lightfield files in directory
            if length(filenames) ~= self.numberFiles
                error = MException( 'LytroDesktop:renameLightfieldFiles:invalidNumberFilenames' ...
                                  , [ 'Number of filenames must be the same of the ' ...
                                    , 'number of ligthfield files in the directory...' ] );
                error.throw();
            end
            
            % Obtain minimum number for files
            orderedFolders          = sort({lytroFiles.folder});
            [~,minimumFolderNumber] = fileparts(orderedFolders{1});
            minimumFolderNumber     = str2double(minimumFolderNumber);
            
            for iFile = 1:self.numberFiles
                % Obtain ordered number of folder
                [rootFolder,folderNumber] = fileparts(lytroFiles(iFile).folder);
                
                % Obtain new folder path
                index = str2double(folderNumber) - minimumFolderNumber + 1;
                newFolderPath = [rootFolder filesep filenames{index}];
                
                % Move file to new folder path
                movefile(lytroFiles(iFile).folder, newFolderPath);
            end
        end
        
        function saveMetadataInSeparateFiles(self)
            %
            % Extract and save metadata from lightfield images in a
            % separate file.
            %
            
            % Obtain lightfield image files metadata
            lytroFiles = self.metadata;

            for iFile = 1:self.numberFiles
                % Define output name for metadata file
                [~,filename] = fileparts(lytroFiles(iFile).name);
                jsonFilepath = [lytroFiles(iFile).folder filesep filename '.json'];
                
                % Write lytro metadata
                LFWriteMetadata(jsonFilepath, lytroFiles(iFile).metadata);
            end
        end
    end
end
