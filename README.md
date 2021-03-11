# US3D
This repository contains somce code to analyze in Matlab 3D ultrasound data in Stradwin format.
The main folder contains two example scripts:

example_stradwin2nifti: example script to convert data in Stradwin format to a NIfTI file. You can open the NIfTI file with, for example, ITK-SNAP (download here: http://www.itksnap.org/pmwiki/pmwiki.php?n=Downloads.SNAP3)
example_sw_to_RCS: example script to read in Stradwin data and build up matrices with room coordinates ('real-world' or 'physical' coordinates) for each pixel of each slice. After running this script, the coordinates are stored in matrices X_RCS, Y_RCS, Z_RCS. The pixel intensities are stored in PXDATA.

The accessory functions are stored in subfolder US3D_functions. The folder misc contains the C++ code used in Stradwin to build up coordinate systems.
