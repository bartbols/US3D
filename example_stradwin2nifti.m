% This script converts 3D ultrasound data in Stradwin format to NIfTI
% format.
clear
% Add required functions to the Matlab path.
addpath(genpath('US3D_functions'))

% Set filename of the stradwin data. It is expected that the corresponding
% .sxi file is stored in the same folder as the .sw file.
filename.sw = 'data/passive.sw';

% Convert the Stradwin data to a NIfTI file.
sw2nifti(filename.sw)