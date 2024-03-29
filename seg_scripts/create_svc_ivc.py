from img import *
import SimpleITK as sitk

import numpy as np
import nrrd
# import pylab
import json
import argparse
import os

parser = argparse.ArgumentParser(description='To run: python3 create_svc_ivc.py [path_to_points]')
parser.add_argument("path_to_points")
args = parser.parse_args()
path2points = args.path_to_points

points_file = open(path2points+'/points.json')
points_data = json.load(points_file)

# ----------------------------------------------------------------------------------------------
# Find the origin and spacing
# ----------------------------------------------------------------------------------------------
# NOTE - We save the origin and spacings because the "save_itk" function used in
# the next step makes the format of the header difficult to work with.
# ----------------------------------------------------------------------------------------------
origin_spacing_file = open(path2points+'/origin_spacing.json')
origin_spacing_data = json.load(origin_spacing_file)
origin = origin_spacing_data['origin']
spacings = origin_spacing_data['spacing']

LV_BP_label = 1
LV_myo_label = 2
RV_BP_label = 3
LA_BP_label = 4
RA_BP_label = 5
Ao_BP_label = 6
PArt_BP_label = 7
LPV1_label = 8
LPV2_label = 9
RPV1_label = 10
RPV2_label = 11
LAA_label = 12
SVC_label = 13
IVC_label = 14

# ----------------------------------------------------------------------------------------------
# Give the path to the segmentation (with pulmonary veins separated)
# ----------------------------------------------------------------------------------------------
seg_corrected_nrrd = path2points+'/seg_corrected.nrrd'

# ----------------------------------------------------------------------------------------------
# Give the paths to the SVC/IVC cylinders and the aorta/pulmonary artery slicers
# Give the associated labels
# ----------------------------------------------------------------------------------------------
svc_nrrd = path2points+'/SVC.nrrd'
ivc_nrrd = path2points+'/IVC.nrrd'

aorta_slicer_nrrd = path2points+'/aorta_slicer.nrrd'
aorta_slicer_label = 0

PArt_slicer_nrrd = path2points+'/PArt_slicer.nrrd'
PArt_slicer_label = 0

# ----------------------------------------------------------------------------------------------
# Convert all of the segmentations to arrays
# ----------------------------------------------------------------------------------------------
seg_s2_array, header1 = nrrd.read(seg_corrected_nrrd)
svc_array, header2 = nrrd.read(svc_nrrd)
ivc_array, header3 = nrrd.read(ivc_nrrd)
aorta_slicer_array, header4 = nrrd.read(aorta_slicer_nrrd)
PArt_slicer_array, header5 = nrrd.read(PArt_slicer_nrrd)

# ----------------------------------------------------------------------------------------------
# Add the SVC and IVC 
# ----------------------------------------------------------------------------------------------
print('\n ## Adding the SVC, IVC and slicers ## \n')
seg_s2a_array = add_masks_replace_only(seg_s2_array, svc_array, SVC_label,RPV1_label)
seg_s2a_array = add_masks(seg_s2a_array, ivc_array, IVC_label)

# ----------------------------------------------------------------------------------------------
# Format and save the segmentation
# ----------------------------------------------------------------------------------------------
print(' ## Formatting and saving the segmentation ##')
seg_s2a_array = np.swapaxes(seg_s2a_array,0,2)
save_itk(seg_s2a_array, origin, spacings, path2points+'/seg_s2a.nrrd')
print(" ## Saved segmentation with SVC/IVC added ##")
