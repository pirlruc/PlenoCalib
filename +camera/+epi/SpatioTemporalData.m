classdef SpatioTemporalData < abstract.TemplateVideo
    %SPATIOTEMPORALDATA
    %   Spatio-temporal data methods and properties.
    
    methods
        function self = SpatioTemporalData(varargin)
            %
            % Spatio-temporal data instance.
            %
            % INPUTS:
            %   1. spatioTemporalData - spatio-temporal data.
            %
            narginchk(0,1);
            
            % Create superclass instance
            self = self@abstract.TemplateVideo();

            if ~isempty(varargin)
                if nargin >= 1
                    self.data = varargin{1};
                end
            end
        end
    end
    
    methods
        function self = addEpipolarPlaneImages_u(self, newEpipolarPlaneImages_u)
            %
            % Add epipolar plane images to spatio-temporal data considering
            % epipolar plane images defined as:
            %           channels x pixel u x frames
            %
            % INPUTS:
            %   1. newEpipolarPlaneImages_u - epipolar plane images to be 
            %   included in the spatio-temporal data.
            %
            narginchk(2,2);
            
            % Add epipolar plane image to spatio-temporal data
            self.data = cat(2,self.data,newEpipolarPlaneImages_u);
        end
        
        function self = addEpipolarPlaneImages_v(self, newEpipolarPlaneImages_v)
            %
            % Add epipolar plane images to spatio-temporal data considering
            % epipolar plane images defined as:
            %           channels x pixel v x frames
            %
            % INPUTS:
            %   1. newEpipolarPlaneImages_v - epipolar plane images to be 
            %   included in the spatio-temporal data.
            %
            narginchk(2,2);
            
            % Add epipolar plane image to epipolar volume
            self.data = cat(3,self.data,newEpipolarPlaneImages_v);
        end
        
        function plot3(self,pixelPosition)
            %
            % Plot spatio-temporal data cutted at a given position.
            %
            % INPUTS:
            %   1. pixelPosition - pixel position to be considered to cut
            %   the spatio-temporal data.
            %
            narginchk(1,2);
            
            % If pixel position is not defined, consider the middle of the
            % spatio-temporal data
            if nargin <= 1
                pixelPosition = [ round(self.size.numberPixels_u / 2) ...
                                ; round(self.size.numberPixels_v / 2) ];
            end
            % Create image pixel instance
            pixelPosition = utils.enums.Classes.PIXEL().convert(pixelPosition);
            
            % Create video instance considering that the frames are given
            % by pixel u:
            %       channels x pivel v x frames x pixel u
            video = camera.Video(permute(self.data,[1,3,4,2]));
            % Obtain image of first frame (x = 1) and cut frame (x = pixel
            % u)
            image_x_1   = image.Image.ImageFromVideo(video,1).data;
            image_x_cut = image.Image.ImageFromVideo(video,pixelPosition.u);

            % Slice the images to obtain the corresponding epipolar plane
            % images orthogonal to the x-axis and cosidering the pixel
            % position for cut.
            image_x_cut = image_x_cut.data(:,1:pixelPosition.v,:);
            
            % Create video instance considering that the frames are given
            % by pixel v:
            %       channels x pivel u x frames x pixel v
            video = camera.Video(permute(self.data,[1,2,4,3]));
            % Obtain image of first frame (y = 1) and cut frame (y = pixel
            % v)
            image_y_1   = image.Image.ImageFromVideo(video,1).data;
            image_y_cut = image.Image.ImageFromVideo(video,pixelPosition.v);

            % Slice the images to obtain the corresponding epipolar plane
            % images orthogonal to the y-axis and cosidering the pixel
            % position for cut.
            image_y_cut = image_y_cut.data(:,1:pixelPosition.u,:);

            % Create video instance using the regular concept for frames:
            %       channels x pixel u x pivel v x frames
            video = camera.Video(self.data);
            % Obtain image of first frame (z = 1) and last frame (z =
            % number frames)
            image_z_1    = image.Image.ImageFromVideo(video,1);
            image_z_last = image.Image.ImageFromVideo(video,video.size.numberFrames).data;
            
            % Slice the image of first frame to consider the pixel position
            % for cut
            image_z_1   = image_z_1.data(:,1:pixelPosition.u,1:pixelPosition.v);
            
            % Encode images to matlab format
            image_x_1    = self.MATLAB_FORMAT.encode(image_x_1);
            image_x_cut  = self.MATLAB_FORMAT.encode(image_x_cut);
            image_y_1    = self.MATLAB_FORMAT.encode(image_y_1);
            image_y_cut  = self.MATLAB_FORMAT.encode(image_y_cut);
            image_z_1    = self.MATLAB_FORMAT.encode(image_z_1);
            image_z_last = self.MATLAB_FORMAT.encode(image_z_last);
            
            % Obtain meshgrids for spatio-temporal data representation
            [image_x_1_x,image_x_1_y,image_x_1_z] = ...
                meshgrid(1,1:self.size.numberPixels_v,1:self.size.numberFrames);
            [image_x_cut_x,image_x_cut_y,image_x_cut_z] = ...
                meshgrid(pixelPosition.u,1:pixelPosition.v,1:self.size.numberFrames);
            [image_y_1_x,image_y_1_y,image_y_1_z] = ...
                meshgrid(1:self.size.numberPixels_u,1,1:self.size.numberFrames);
            [image_y_cut_x,image_y_cut_y,image_y_cut_z] = ...
                meshgrid(1:pixelPosition.u,pixelPosition.v,1:self.size.numberFrames);
            [image_z_1_x,image_z_1_y,image_z_1_z] = ...
                meshgrid(1:pixelPosition.u,1:pixelPosition.v,1);
            [image_z_last_x,image_z_last_y,image_z_last_z] = ...
                meshgrid(1:self.size.numberPixels_u,1:self.size.numberPixels_v,self.size.numberFrames);
            
            % The meshgrid function in matlab defines arrays with the
            % internal format switched in the first two dimensions. This
            % means, that to get the internal format we shoud switch the
            % columns by the rows. Since we want the meshgrids in matlab
            % format, we do not change the columns by the rows. We only
            % consider this switch to eliminate the dimension that only has
            % one item.
            image_x_1_x   = permute(image_x_1_x,[3,1,2]);
            image_x_1_y   = permute(image_x_1_y,[3,1,2]);
            image_x_1_z   = permute(image_x_1_z,[3,1,2]);
            image_x_cut_x = permute(image_x_cut_x,[3,1,2]);
            image_x_cut_y = permute(image_x_cut_y,[3,1,2]);
            image_x_cut_z = permute(image_x_cut_z,[3,1,2]);
            image_y_1_x   = permute(image_y_1_x,[3,2,1]);
            image_y_1_y   = permute(image_y_1_y,[3,2,1]);
            image_y_1_z   = permute(image_y_1_z,[3,2,1]);
            image_y_cut_x = permute(image_y_cut_x,[3,2,1]);
            image_y_cut_y = permute(image_y_cut_y,[3,2,1]);
            image_y_cut_z = permute(image_y_cut_z,[3,2,1]);
            
            % Plot volume using slices at x = 1, y = 1, z = number frames
            surface( image_x_1_x,image_x_1_y,image_x_1_z ...
                   , image_x_1 ...
                   , 'FaceColor', 'texturemap' ...
                   , 'EdgeColor', 'none' ...
                   , 'CDataMapping', 'direct' );
            surface( image_y_1_x,image_y_1_y,image_y_1_z ...
                   , image_y_1 ...
                   , 'FaceColor', 'texturemap' ...
                   , 'EdgeColor', 'none' ...
                   , 'CDataMapping', 'direct' );
            surface( image_z_last_x,image_z_last_y,image_z_last_z ...
                   , image_z_last ...
                   , 'FaceColor', 'texturemap' ...
                   , 'EdgeColor', 'none' ...
                   , 'CDataMapping', 'direct' );
               
            % Plot volume slice using slices at x = pixel u, y = pixel v
            % and z = 1
            surface( image_x_cut_x,image_x_cut_y,image_x_cut_z ...
                   , image_x_cut ...
                   , 'FaceColor', 'texturemap' ...
                   , 'EdgeColor', 'none' ...
                   , 'CDataMapping', 'direct' );
            surface( image_y_cut_x,image_y_cut_y,image_y_cut_z ...
                   , image_y_cut ...
                   , 'FaceColor', 'texturemap' ...
                   , 'EdgeColor', 'none' ...
                   , 'CDataMapping', 'direct' );
            surface( image_z_1_x,image_z_1_y,image_z_1_z ...
                   , image_z_1 ...
                   , 'FaceColor', 'texturemap' ...
                   , 'EdgeColor', 'none' ...
                   , 'CDataMapping', 'direct' );
            xlabel('Pixel u');
            ylabel('Pixel v');
            zlabel('Frames');
        end
    end
    
    methods (Static)
        function self = SpatioTemporalDataFromFile(videoFilepath)
            %
            % Create spatio-temporal data instance from file.
            %
            % INPUTS:
            %   1. videoFilepath - video local filepath.
            %
            narginchk(1,1);
            
            % Read video data. Convert data to double.
            videoData = im2double(VideoReader(videoFilepath).read);
            
            % Obtain matlab format and decode
            matlabFormat = camera.enums.Formats.MATLAB_FORMAT();
            videoData    = matlabFormat.decode(videoData);
            
            % Create video instance
            self = camera.epi.SpatioTemporalData(videoData);
        end
    end
end
