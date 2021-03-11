% This script converts 3D ultrasound data in Stradwin format to NIfTI
% format.
clear
% Add required functions to the Matlab path.
addpath(genpath('US3D_functions'))

% Set filename of the stradwin data. It is expected that the corresponding
% .sxi file is stored in the same folder as the .sw file.
filename.sw = 'data/passive.sw';

% Set to true if you want to plot the images and landmarks in 3D. By
% default, only every 5th slice is shown.
PlotFlag = true;
%% Build up calibration transformation matrix.
%% Read the Stradwin data.
[PXDATA,RES,IM,LANDMARK,OBJECT] = read_sw(filename.sw);
% Get x- and y-scale (pixel dimensions)
xscale = str2double(RES.RES_XSCALE); % cm per px
yscale = str2double(RES.RES_YSCALE); % cm per px

% See Stradwin website for more information:
% http://mi.eng.cam.ac.uk/~gmt11/stradwin/stradwin_files.htm
T_cal = make_transform_matrix(...
    str2double(RES.RES_XTRANS),...
    str2double(RES.RES_YTRANS),...
    str2double(RES.RES_ZTRANS),...
    str2double(RES.RES_AZIMUTH),...
    str2double(RES.RES_ELEVATION),...
    str2double(RES.RES_ROLL),'inverse');
%% Transform to room coordinate system
% Set pixel_skip and slice_skip to values different from 1 to read in
% downsampled data using only the n-th pixel and slice (e.g. slice_skip = 2
% reads in every other slice, pixel_skip = 2 reads in 2D US images at half 
% the resolution of the original image).

pixel_skip = 1;
slice_skip = 1;

% Create a matrix with pixel coordinates of all pixels in the image.
% The 0.5 is used because the center of the corner pixel is half a pixel 
% away from the edge of the image.
w = size(PXDATA,2);
h = size(PXDATA,1);


[Xp,Yp] = meshgrid(0.5:pixel_skip:w-0.5,0.5:pixel_skip:h-0.5); 
Xp_vec = Xp(:) * xscale;
Yp_vec = Yp(:) * yscale;

PXDATA = PXDATA(1:pixel_skip:end,1:pixel_skip:end,1:slice_skip:end);

% Coordinates of image corners.
Xc = [0 0 w*xscale w*xscale]';
Yc = [0 h*yscale 0 h*yscale]';

if PlotFlag == true
    figure('Color','w')
    hold on
    clim = [0 max(PXDATA(:))];
end

% Create empty matrices for all pixel coordinates.
X_RCS = zeros(size(PXDATA));
Y_RCS = zeros(size(PXDATA));
Z_RCS = zeros(size(PXDATA));
% AllPX = zeros(numel(PXDATA),1);
first = 1;
for frame_nr = 1 : 1 : size(PXDATA,3)
    % Transform from image to room coordinates, using the description
    % of the Stradwin coordinate system:
    % http://mi.eng.cam.ac.uk/~gmt11/stradwin/stradwin_files.htm
    T_frame = make_transform_matrix(...
        IM(frame_nr).values(1),...
        IM(frame_nr).values(2),...
        IM(frame_nr).values(3),...
        IM(frame_nr).values(4),...
        IM(frame_nr).values(5),...
        IM(frame_nr).values(6),'inverse');
    
    T_iso = make_transform_matrix(...
        str2double(RES.RES_ISOCENTRE_XTRANS),...
        str2double(RES.RES_ISOCENTRE_YTRANS),...
        str2double(RES.RES_ISOCENTRE_ZTRANS),...
        str2double(RES.RES_ISOCENTRE_AZIMUTH),...
        str2double(RES.RES_ISOCENTRE_ELEVATION),...
        str2double(RES.RES_ISOCENTRE_ROLL),'inverse');
    
    XYZp_room = (T_cal * T_frame * T_iso) \ [Xp_vec Yp_vec zeros(size(Xp_vec(:))) ones(size(Xp_vec(:)))]';
    XYZc_room = (T_cal * T_frame * T_iso) \ [Xc Yc zeros(size(Xc)) ones(size(Xc))]';
    
    % Get pixel data.
    frame_data = squeeze(PXDATA(:,:,frame_nr));
    %     hs = scatter3(XYZp_room(1,:),XYZp_room(2,:),XYZp_room(3,:),3,frame_data(:)');
    
    % Store all data in separate vectors to usedlater with the gridded 
    % interpolant function.
    last = first + w*h - 1;
    X_RCS(:,:,frame_nr)  = reshape(XYZp_room(1,:),size(frame_data));
    Y_RCS(:,:,frame_nr)  = reshape(XYZp_room(2,:),size(frame_data));
    Z_RCS(:,:,frame_nr)  = reshape(XYZp_room(3,:),size(frame_data));
%     AllPX(first:last,1) = frame_data(:);
    first = last + 1;
    
    if PlotFlag == true && mod(frame_nr,5)==0 % show every 5th slice
      % Create plot objects
        himg = image([0 w*xscale],[0 h*yscale],frame_data);
        hp = patch('Vertices',[Xc Yc zeros(size(Xc))],...
            'Faces',[1 2 4 3],...
            'EdgeColor','none',...
            'FaceColor','none',...
            'LineWidth',1);
        
        % Set transforms on the plot objects to display the images in the Stradwin
        % room coordinate system
        t = hgtransform('Matrix',inv(T_cal * T_frame * T_iso));
        set([hp,himg],'Parent',t)
%         set([hp,himg],'Parent',t)
        colormap gray
        axis equal tight on
        set(gca,'CLim',clim);
    end
end

% Plot the landmarks
if PlotFlag == true
    for nr = 1 : length(LANDMARK)
        plot3(LANDMARK(nr).values(1),...
            LANDMARK(nr).values(2),...
            LANDMARK(nr).values(3),...
            'o','MarkerFaceColor','r',...
            'MarkerEdgeColor','k',...
            'MarkerSize',5)
        text(LANDMARK(nr).values(1),...
            LANDMARK(nr).values(2),...
            LANDMARK(nr).values(3),...
            int2str(LANDMARK(nr).nr))
        
    end
    xlabel('x');ylabel('y');zlabel('z')
    set(gca,'FontSize',14)
end
