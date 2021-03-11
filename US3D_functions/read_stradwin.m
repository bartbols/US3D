clear
filename.sw = 'data/passive.sw';
PlotFlag = true;
[PXDATA,RES,IM,LANDMARK,OBJECT] = read_sw(filename.sw);
% f = RES_BUF_FRAMES 287
% w = RES_BUF_WIDTH 513
% h = RES_BUF_HEIGHT 480

xscale = str2double(RES.RES_XSCALE); % cm per px
yscale = str2double(RES.RES_YSCALE); % cm per px

T_cal = make_transform_matrix(...
    str2double(RES.RES_XTRANS),...
    str2double(RES.RES_YTRANS),...
    str2double(RES.RES_ZTRANS),...
    str2double(RES.RES_AZIMUTH),...
    str2double(RES.RES_ELEVATION),...
    str2double(RES.RES_ROLL),'inverse');
%%
% Create a matrix with pixel coordinates of all pixels in the image.
% The 0.5 is used because the center of the corner pixel is half a pixel 
% away from the edge of the image.
w = size(PXDATA,2);
h = size(PXDATA,1);
pixel_skip = 1;
slice_skip = 1;
[Xp,Yp] = meshgrid(0.5:pixel_skip:w-0.5,0.5:pixel_skip:h-0.5); 
Xp_vec = Xp(:) * xscale;
Yp_vec = Yp(:) * yscale;
PXDATA = PXDATA(1:pixel_skip:end,1:pixel_skip:end,1:slice_skip:end);

% PXDATA = PXDATA(:,:,100:150);
% IM = IM(100:150);
% Coordinates of image corners.
Xc = [0 0 w*xscale w*xscale]';
Yc = [0 h*yscale 0 h*yscale]';

if PlotFlag == true
    figure('Color','w')
    hold on
    clim = [0 max(PXDATA(:))];
end

% Create empty matrices for all pixel coordinates.
AllX = zeros(size(PXDATA));
AllY = zeros(size(PXDATA));
AllZ = zeros(size(PXDATA));
% AllPX = zeros(numel(PXDATA),1);
first = 1;
for frame_nr = 1 : 5 : size(PXDATA,3)
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
    
    % Store all data in separate vectors to used later with the gridded 
    % interpolant function.
    last = first + w*h - 1;
    AllX(:,:,frame_nr)  = reshape(XYZp_room(1,:),size(frame_data));
    AllY(:,:,frame_nr)  = reshape(XYZp_room(2,:),size(frame_data));
    AllZ(:,:,frame_nr)  = reshape(XYZp_room(3,:),size(frame_data));
%     AllPX(first:last,1) = frame_data(:);
    first = last + 1;
    
    if PlotFlag == true
        % Create plot objects
        himg = image([0 w*xscale],[0 h*yscale],frame_data);
        hp = patch('Vertices',[Xc Yc zeros(size(Xc))],...
            'Faces',[1 2 4 3],...
            'EdgeColor','y',...
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

%% Interpolate to a regular grid.
tic
F = scatteredInterpolant(AllX(:),AllY(:),AllZ(:),double(PXDATA(:)),...
    'linear','none');
t_elapsed = toc;

fprintf('It took %.2f seconds to build up the interpolant function.\n',t_elapsed)

%% Interpolate to a regular 3D grid. 
tic
% x_range = min(AllX(:)) : gridsize(1) : max(AllX(:));
% y_range = min(AllY(:)) : gridsize(2) : max(AllY(:));
% z_range = min(AllZ(:)) : gridsize(3) : max(AllZ(:));
gridsize = [0.02 0.02 0.02]; % in cm
x_range = [ 12 : gridsize(1) : 14];
y_range = [-19 : gridsize(2) : -17];
z_range = [-29 : gridsize(3) : -27];


[Xg,Yg,Zg] = ndgrid(...
    x_range,...
    y_range,...
    z_range);

IMG_3D = F(Xg,Yg,Zg);
IMG_3D(isnan(IMG_3D)) = 0;
t_elapsed = toc;

fprintf('It took %.2f seconds to interpolate a 3D image.\n',t_elapsed)
% Write as a NIfTI file
load info
info.PixelDimensions = gridsize;
info.ImageSize = size(IMG_3D);
% niftiwrite(IMG_3D,'test.nii',info,'Compressed',true)
niftiwrite(IMG_3D,'passive.nii',info,'Compressed',true)

%% Apply 3D difference of Gaussian filter
sigma1 = 4;
% sigma2 = 5;
filter_size = 21;
h1 = fspecial3('gaussian',filter_size,sigma1);
% h2 = fspecial3('gaussian',filter_size,sigma2);
diff_filter = h1-h2;

gLeft = zeros(filter_size,filter_size+1,filter_size);
gLeft(:,1:filter_size,:)  = h1;
gRight = zeros(filter_size,filter_size+1,filter_size);
gRight(:,2:filter_size+1,:)  = h1;
gDiff = gLeft-gRight;
DoG.Gx = gDiff(:,1:filter_size,:);

gTop = zeros(filter_size+1,filter_size,filter_size);
gTop(1:filter_size,:,:)  = h1;
gBottom = zeros(filter_size+1,filter_size,filter_size);
gBottom(2:filter_size+1,:,:)  = h1;
gDiff = gTop-gBottom;
DoG.Gy = gDiff(1:filter_size,:,:);

gUp = zeros(filter_size,filter_size,filter_size+1);
gUp(:,:,1:filter_size)  = h1;
gDown = zeros(filter_size,filter_size,filter_size+1);
gDown(:,:,2:filter_size+1)  = h1;
gDiff = gUp-gDown;
DoG.Gz = gDiff(:,:,1:filter_size);

% gDiff = gLeft-gRight;
% DoG.Gx = gDiff(:,1:maskSize);
Ix = convn(IMG_3D,DoG.Gx,'same');
Iy = convn(IMG_3D,DoG.Gy,'same');
Iz = convn(IMG_3D,DoG.Gz,'same');

% niftiwrite(im_filt,'filt.nii',info,'Compressed',true)

% %%
% % [rows, cols,slices] = size(IMG_3D);
% % maskSize = max([rows, cols,slices]); 
% % midpt = ceil(maskSize/2);
% DoG = difference_of_gaussian_kernels(11);
% % 
% Ix = convn(IMG_3D, DoG.Gx,'same');
% Iy = convn(IMG_3D, DoG.Gy,'same');
% Iz = convn(IMG_3D, DoG.Gy,'same');
niftiwrite(Ix,'Ix.nii',info,'Compressed',true)
niftiwrite(Iy,'Iy.nii',info,'Compressed',true)
niftiwrite(Iz,'Iz.nii',info,'Compressed',true)

%% Make test data
IMG_3D = zeros(51,51,51);
IMG_3D(1:3:end,1:3:end,:) = 1;
% IMG_3D(1:3:end,:,1:3:end) = 1;
% IMG_3D = IMG_3D + randn(size(IMG_3D))*0.1;
% [Gx,Gy,Gz] = imgradientxyz(IMG_3D);
sigma = 1;
% sigma1 = 2;
% sigma2 = 3;

% im_filt = imgaussfilt3(IMG_3D,sigma);

% Difference of gaussian image. Not sure if this is correct...
% im_filt = imgaussfilt3(IMG_3D,sigma1) - imgaussfilt3(IMG_3D,sigma2) ;

[Ix,Iy,Iz] = imgradientxyz(IMG_3D);
% [Ix,Iy,Iz] = imgradientxyz(im_filt);

load info
info.PixelDimensions = gridsize;
info.ImageSize = size(im_filt);
niftiwrite(im_filt,'filt.nii',info,'Compressed',true)
%% Compare image gradient with partial_derivative
[Ix2,Iy2,Iz2] = partial_derivative_3D(IMG_3D);
Sdata = partial_derivative_to_structure_tensor_form([Ix2 Iy2 Iz2]);
[e1,e2,e3,l1,l2,l3] = eigen_decomposition(Sdata);
%% Write NIfTI file with gradient direction
[Gmag,Gazimuth,Gelevation] = imgradient3(Ix,Iy,Iz);
% GradDir = 
[x,y,z] = sph2cart(Gazimuth,Gelevation,1);
GradDir = zeros([size(IMG_3D) 1 3]);
GradDir(:,:,:,1,1) = x;
GradDir(:,:,:,1,2) = y;
GradDir(:,:,:,1,3) = z;
load info
info.ImageSize = size(GradDir);
info.PixelDimensions = [gridsize 0 0];
niftiwrite(GradDir,'graddir.nii',info,'Compressed',true)

%% Calculate structure tensor
E1 = zeros([size(IMG_3D) 3]);
E2 = zeros([size(IMG_3D) 3]);
E3 = zeros([size(IMG_3D) 3]);
L1 = zeros(size(IMG_3D));
L2 = zeros(size(IMG_3D));
L3 = zeros(size(IMG_3D));
for i = 1 : size(Ix,1)
    tic
    for j = 1 : size(Ix,2)
        for k = 1 : size(Ix,3)
%             if IMG_3D(i,j,k) == 0;continue;end

            S = [Ix(i,j,k)^2         Ix(i,j,k)*Iy(i,j,k) Ix(i,j,k)*Iz(i,j,k);...
                 Ix(i,j,k)*Iy(i,j,k) Iy(i,j,k)^2         Iy(i,j,k)*Iz(i,j,k);...
                 Ix(i,j,k)*Iz(i,j,k) Iy(i,j,k)*Iz(i,j,k) Iz(i,j,k)^2        ];
             
              [e1,e2,e3,l1,l2,l3] = eigen_decomposition(S);
             
              E1(i,j,k,1:3) = e1;
              E2(i,j,k,1:3) = e2;
              E3(i,j,k,1:3) = e3;
              L1(i,j,k) = l1;
              L2(i,j,k) = l2;
              L3(i,j,k) = l3;
        end
    end
    toc
end
FA = FuncFA( L1,L2,L3 );
load info
info.ImageSize = size(E1);
info.PixelDimensions = [gridsize 0];
niftiwrite(abs(E1),'E1.nii',info,'Compressed',true)
niftiwrite(abs(E2),'E2.nii',info,'Compressed',true)
niftiwrite(abs(E3),'E3.nii',info,'Compressed',true)

info.ImageSize = size(FA);
info.PixelDimensions = gridsize;
niftiwrite(FA,'FA.nii',info,'Compressed',true)

%%
k = 170;
figure('Color','w')
subplot(2,2,1);hold on
% image(IMG_3D(:,:,k))
image(im_filt(:,:,k))
colormap gray
subplot(2,2,2);hold on
image(squeeze(E1(:,:,k,:)))
% imagesc(FA(:,:,k))
% imagesc(L1(:,:,k))
subplot(2,2,3);hold on
imagesc(L2(:,:,k))
subplot(2,2,4);hold on
imagesc(L3(:,:,k))
linkprop(findobj(gcf,'Type','Axes'),{'XLim','YLim'})
