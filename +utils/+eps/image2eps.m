function image2eps( currentFolder )
%IMAGE2EPS
%   Convert image files to eps files.

    if nargin < 1
        currentFolder = pwd;
    end 
    
    % Obtain filenames 
    filenames = dir(currentFolder);
    filenames = filenames(~[filenames.isdir]);

    for iFilename = 1:length(filenames)
        % Obtain filepath
        filepath = [filenames(iFilename).folder filesep filenames(iFilename).name];
        
        % Create image instance
        try
            oImage = image.Image.ImageFromFile(filepath);
        catch
            continue
        end
        
        % Obtain new filename with eps extension
        [path,filename] = fileparts(filepath);
        outputFilepath  = [path filesep filename '.eps'];
        
        % Write image file to eps file format
        oImage.write(outputFilepath);
    end
end
