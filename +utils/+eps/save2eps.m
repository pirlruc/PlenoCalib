function save2eps( outputFilename, options )
%SAVE2EPS
%   Save figure to eps file.
    
    if nargin < 1
        outputFilename = 'outputImage.eps';
    end
    
    if nargin < 2
        options = [];
    end

    % Set default text fonts and sizes to save file
    REMOVE_TITLE = true;
    FONT_SIZE    = 24;
    utils.vislab.toprint(gcf,REMOVE_TITLE,FONT_SIZE,options);
    
    % Set axis options to remove extra white space in figure
    set(gca,'LooseInset',get(gca,'TightInset')); 
    
    % Save eps color file
    saveas(gcf,outputFilename,'epsc');
end
