function sw2nifti(filename_sw,varargin)
%SW2NIFTI Reads in 3D ultrasound data in NIfTI format and writes the data
%as a NIfTI file.
%
% Bart Bolsterlee
% Neuroscience Research Australia
% March 2021
%
% ----------------- USAGE -----------------
% sw2nifti( filename);
% sw2nifti( filename,'filename_nii','test.nii');
% 
% ----------------- INPUT -----------------
% - filename_sw       :  Filename of a Stradwin file containing the
%                        metadata for the ultrasound images. It is expected
%                        that the corresponding .sxi file lives in the same
%                        folder as the .sw file.
%
% Optional input, provided as 'parameter',<value> pairs:
% - filename_nii : filename of the NIfTI file. Default: same name as .sw
%                  file but with extension .nii (or .nii.gz if 'compressed'=true)
% - compressed   : logical flag (true/false) indicating if NIfTI file is
%                  compressed (.nii.gz) or not. Default: true
% - slice_spacing: spacing in mm between ultrasound images. Default: 1
%         
% ----------------- OUTPUT -----------------
% this function does not have outputs, but a NIfTI file will be created.


% Read inputs to function
p = inputParser;
addRequired(p,'filename_sw')
addParameter(p,'filename_nii',[],@(x) contains('.nii'))
addParameter(p,'compressed',true,@(x) x==0 || x==1 || islogical(x) )
addParameter(p,'slice_spacing',1,@(x) isscalar(x) && x>=0)
parse(p, filename_sw,varargin{:})

if isempty(p.Results.filename_nii)
    filename_nii = strrep(filename_sw,'.sw','.nii');
else
    filename_nii = p.Results.filename_nii;
end

%% Read the stradwin data
fprintf('Reading data from %s... ',filename_sw)
[PXDATA,RES] = read_sw(filename_sw);
fprintf(' completed.\n')

% Get x- and y-scale (pixel dimensions)
xscale = str2double(RES.RES_XSCALE); % cm per px
yscale = str2double(RES.RES_YSCALE); % cm per px

% Write as NIfTI file
% load info_int16
% Load NIfTI header structure.
load(fullfile(fileparts(mfilename('fullpath')),'info_int16.mat'))
info.PixelDimensions = [xscale*10 yscale*10 p.Results.slice_spacing];
info.ImageSize = size(PXDATA);
fprintf('Writing %d ultrasound slices to %s',size(PXDATA,3),filename_nii)
if p.Results.compressed
    fprintf('.gz')
end
niftiwrite(PXDATA,strrep(filename_nii,'.sw','.nii'),...
    info,'Compressed',p.Results.compressed)
fprintf('. Completed.\n')

end

