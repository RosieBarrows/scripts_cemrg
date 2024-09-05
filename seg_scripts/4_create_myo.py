import cut_labels
from img import *

import SimpleITK as sitk
import numpy as np
import nrrd
import json
import os
import argparse

#----------------------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description='To run: python3 create_myo.py [path_to_points]')
parser.add_argument("path_to_points")
args = parser.parse_args()
path2points = args.path_to_points

if os.path.isfile(f"{path2points}/points.txt"):
	os.system("python3 txt_2_json.py "+path2points+"/points.txt "+path2points+"/labels.txt "+path2points+"/points.json")

points_file = open(path2points+'/points.json')
points_data = json.load(points_file)

thickness_file = open(f"{path2points}/thickness.json")
thickness_data = json.load(thickness_file)

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

# ----------------------------------------------------------------------------------------------
# Make a "tmp" folder for temporarily storing steps in the segmentation pipeline
# ----------------------------------------------------------------------------------------------
os.system("mkdir -p "+path2points+"/tmp")
# ----------------------------------------------------------------------------------------------
# Define the wall thickness
# ----------------------------------------------------------------------------------------------
# sf = 1/0.39844 # scale factor
sf = np.ceil(1.0/spacings[0]) # scale factor

LV_neck_WT = sf*thickness_data["aorta"]
RV_WT = sf*thickness_data["RV"]
LA_WT = sf*thickness_data["LA"]
RA_WT = sf*thickness_data["RA"]
Ao_WT = sf*thickness_data["aorta"]
PArt_WT = sf*thickness_data["RV"]
rings_thickness = sf*thickness_data["rings"]

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
LV_neck_label = 101
RV_myo_label = 103
LA_myo_label = 104
RA_myo_label = 105
Ao_wall_label = 106
PArt_wall_label = 107

seg_array_good_header = sitk.ReadImage(path2points+'/seg_s2a.nrrd')

# ----------------------------------------------------------------------------------------------
# Prepare the seeds for the tips of the aorta and pulmonary artery
# ----------------------------------------------------------------------------------------------
# Ao_tip_seed = points_data['Ao_tip']
# PArt_tip_seed = points_data['PArt_tip']

# ----------------------------------------------------------------------------------------------
# Connected component in the aorta and save the segmentation
# ----------------------------------------------------------------------------------------------
# print(' ## Running connected component in the aorta ## \n')
# seg_s2f_nrrd = path2points+'/seg_s2f.nrrd'
# seg_aorta_cc_array = connected_component(seg_s2f_nrrd, Ao_tip_seed, Ao_BP_label,path2points)
# seg_aorta_cc_array = np.swapaxes(seg_aorta_cc_array,0,2)
# save_itk_keeping_header(new_image=seg_aorta_cc_array, original_image=seg_array_good_header, filename=path2points+'/seg_aorta_cc.nrrd')


# ----------------------------------------------------------------------------------------------
# Connected component in the pulmonary artery and save the segmentation
# ----------------------------------------------------------------------------------------------
# print(' ## Running connected component in the pulmonary artery ## \n')
# seg_aorta_cc_nrrd = path2points+'/seg_aorta_cc.nrrd'
# seg_PA_cc_array = connected_component(seg_aorta_cc_nrrd, PArt_tip_seed, PArt_BP_label,path2points)
# seg_PA_cc_array = np.swapaxes(seg_PA_cc_array,0,2)

# save_itk_keeping_header(new_image=seg_PA_cc_array, original_image=seg_array_good_header, filename=path2points+'/seg_PA_cc.nrrd')


# ----------------------------------------------------------------------------------------------
# Create the LV outflow tract myocardium
# ----------------------------------------------------------------------------------------------
print('\n ## Step 1/10: Creating myocardium for the LV outflow tract ## \n')
print(' ## LV neck: Executing distance map ## \n')
LV_DistMap = distance_map(path2points+'/seg_s2f.nrrd',LV_BP_label)
print(' ## LV neck: Writing temporary image ## \n')
sitk.WriteImage(LV_DistMap,path2points+'/tmp/LV_DistMap.nrrd',True)

print(' ## LV neck: Thresholding distance filter ## \n')
LV_neck = threshold_filter_nrrd(path2points+'/tmp/LV_DistMap.nrrd',0,LV_neck_WT)
sitk.WriteImage(LV_neck,path2points+'/tmp/LV_neck.nrrd',True)

print(' ## LV neck: Adding LV neck to LV myo ## \n')
LV_neck_array, header = nrrd.read(path2points+'/tmp/LV_neck.nrrd')
seg_s2f_array, header = nrrd.read(path2points+'seg_s2f.nrrd')
LV_neck_array = add_masks_replace(LV_neck_array,LV_neck_array,LV_neck_label)
seg_s3a_array = add_masks(seg_s2f_array,LV_neck_array,2)

# ----------------------------------------------------------------------------------------------
# Format and save the segmentation
# ----------------------------------------------------------------------------------------------
print(' ## LV neck: Formatting and saving the segmentation ## \n')
seg_s3a_array = np.swapaxes(seg_s3a_array,0,2)
# save_itk(seg_s3a_array, origin, spacings, path2points+'/seg_s3a.nrrd')
save_itk_keeping_header(new_image=seg_s3a_array, original_image=seg_array_good_header, filename=path2points+'/seg_s3a.nrrd')
print(" ## LV neck: Saved segmentation with LV outflow tract added ## \n")

# ----------------------------------------------------------------------------------------------
# LV_myo is pushed by RV_BP
# ----------------------------------------------------------------------------------------------
print(' ## Pushing LV_myo with RV_BP ## \n')
seg_s3a_array = push_inside(path2points,path2points+'seg_s3a.nrrd',RV_BP_label,LV_myo_label,LV_BP_label,LV_neck_WT)
seg_s3a_array = np.swapaxes(seg_s3a_array,0,2)
# save_itk(seg_s3a_array, origin, spacings, path2points+'/seg_s3a.nrrd')
save_itk_keeping_header(new_image=seg_s3a_array, original_image=seg_array_good_header, filename=path2points+'/seg_s3a.nrrd')
print(" ## LV neck: Preventing hole in neck due to RV BP ## \n")

# ----------------------------------------------------------------------------------------------
# Create the aortic wall
# ----------------------------------------------------------------------------------------------
print(' ## Step 2/10: Creating the aortic wall ## \n')
print(' ## Aortic wall: Executing distance map ## \n')
Ao_DistMap = distance_map(path2points+'seg_s3a.nrrd',Ao_BP_label)
print(' ## Aortic wall: Writing temporary image ## \n')
sitk.WriteImage(Ao_DistMap,path2points+'/tmp/Ao_DistMap.nrrd',True)

print(' ## Aortic wall: Thresholding distance filter ## \n')
Ao_wall = threshold_filter_nrrd(path2points+'/tmp/Ao_DistMap.nrrd',0,Ao_WT)
sitk.WriteImage(Ao_wall,path2points+'/tmp/Ao_wall.nrrd',True)

print(' ## Aortic wall: Adding aortic wall to segmentation ## \n')
Ao_wall_array, header = nrrd.read(path2points+'/tmp/Ao_wall.nrrd')
seg_s3a_array, header = nrrd.read(path2points+'seg_s3a.nrrd')
Ao_wall_array = add_masks_replace(Ao_wall_array,Ao_wall_array,Ao_wall_label)
seg_s3b_array = add_masks_replace_except_2(seg_s3a_array,Ao_wall_array,Ao_wall_label,LV_BP_label,LV_myo_label)

# ----------------------------------------------------------------------------------------------
# Format and save the segmentation
# ----------------------------------------------------------------------------------------------
print(' ## Aortic wall: Formatting and saving the segmentation ## \n')
seg_s3b_array = np.swapaxes(seg_s3b_array,0,2)
# save_itk(seg_s3b_array, origin, spacings, path2points+'/seg_s3b.nrrd')
save_itk_keeping_header(new_image=seg_s3b_array, original_image=seg_array_good_header, filename=path2points+'/seg_s3b.nrrd')
print(" ## Aortic wall: Saved segmentation with aortic wall added ## \n")

print("## Opening aorta ##\n")

cut_labels.open_artery(path2image      = f"{path2points}/seg_s3b.nrrd", 
					   myo_label       = Ao_wall_label, 
					   artery_bp_label = Ao_BP_label, 
					   ventricle_label = LV_BP_label, 
					   cut_ratio       = 0.95,
					   save_filename   = f"{path2points}/seg_s3b_aorta.nrrd")

# cut_labels.open_artery(path2image      = f"{path2points}/seg_s3d_aorta.nrrd", 
# 					   myo_label       = PArt_wall_label, 
# 					   artery_bp_label = PArt_BP_label, 
# 					   ventricle_label = RV_BP_label, 
# 					   cut_ratio       = 0.85,
# 					   save_filename   = f"{path2points}/seg_s3e.nrrd")


# ----------------------------------------------------------------------------------------------
# Create the pulmonary artery wall
# ----------------------------------------------------------------------------------------------
print(' ## Step 3/10: Creating the pulmonary artery wall ## \n')
print(' ## Pulmonary artery wall: Executing distance map ## \n')
PArt_DistMap = distance_map(f"{path2points}/seg_s3b_aorta.nrrd",PArt_BP_label)
print(' ## Pulmonary artery wall: Writing temporary image ## \n')
sitk.WriteImage(PArt_DistMap,path2points+'/tmp/PArt_DistMap.nrrd',True)

print(' ## Pulmonary artery wall: Thresholding distance filter ## \n')
PArt_wall = threshold_filter_nrrd(path2points+'/tmp/PArt_DistMap.nrrd',0,PArt_WT)
sitk.WriteImage(PArt_wall,path2points+'/tmp/PArt_wall.nrrd',True)

print(' ## Pulmonary artery wall: Adding pulmonary artery wall to segmentation ## \n')
PArt_wall_array, header = nrrd.read(path2points+'/tmp/PArt_wall.nrrd')
seg_s3b_array, header = nrrd.read(f"{path2points}/seg_s3b_aorta.nrrd")
PArt_wall_array = add_masks_replace(PArt_wall_array,PArt_wall_array,PArt_wall_label)
# The pulmonary artery doesn't modify the RV bloodpool or the wall of the aorta, but can modify the aorta.
seg_s3c_array = add_masks_replace_except_3(imga = seg_s3b_array,
										 imgb = PArt_wall_array,
										 newmask = PArt_wall_label,
										 forbid_change1 = RV_BP_label,
										 forbid_change2 = Ao_wall_label,
										 forbid_change3 = Ao_BP_label)
# seg_s3c_array = add_masks_replace_except_2(seg_s3b_array,PArt_wall_array,PArt_wall_label,3,Ao_wall_label)

# ----------------------------------------------------------------------------------------------
# Format and save the segmentation
# ----------------------------------------------------------------------------------------------
print(' ## Pulmonary artery wall: Formatting and saving the segmentation ## \n')
seg_s3c_array = np.swapaxes(seg_s3c_array,0,2)
# save_itk(seg_s3c_array, origin, spacings, path2points+'/seg_s3c.nrrd')
save_itk_keeping_header(new_image=seg_s3c_array, original_image=seg_array_good_header, filename=path2points+'/seg_s3c.nrrd')

print(" ## Pulmonary artery wall: Saved segmentation with pulmonary artery wall added ## \n")

print("## Opening PA ##\n")

cut_labels.open_artery(path2image      = f"{path2points}/seg_s3c.nrrd", 
					   myo_label       = PArt_wall_label, 
					   artery_bp_label = PArt_BP_label, 
					   ventricle_label = RV_BP_label, 
					   cut_ratio       = 0.85,
					   save_filename   = f"{path2points}/seg_s3c_PA.nrrd")
# ----------------------------------------------------------------------------------------------
# Push the pulmonary artery wall into the pulmonary artery blood pool
# ----------------------------------------------------------------------------------------------
# print(' ## Pulmonary artery wall: Pushing the wall of the pulmonary artery ## \n')
# seg_s3d_array = push_inside(path2points,path2points+'seg_s3c.nrrd',Ao_wall_label,PArt_wall_label,PArt_BP_label,PArt_WT)

# print(' ## Pulmonary artery wall: aorta pushes PA ## \n')
# seg_s3d_array = push_inside(path2points = path2points,
# 							img_nrrd = f"{path2points}/seg_s3c_PA.nrrd",pusher_wall_lab = PArt_wall_label,
# 							pushed_wall_lab = Ao_wall_label,
# 							pushed_BP_lab = Ao_BP_label,
# 							pushed_WT = Ao_WT)

# seg_s3d_array = push_inside(path2points = path2points,
# 							img_nrrd = f"{path2points}/seg_s3c_PA.nrrd",pusher_wall_lab = Ao_wall_label,
# 							pushed_wall_lab = PArt_wall_label,
# 							pushed_BP_lab = PArt_BP_label,
# 							pushed_WT = PArt_WT)

# ----------------------------------------------------------------------------------------------
# Format and save the segmentation
# ----------------------------------------------------------------------------------------------
# print(' ## Pulmonary artery wall: Formatting and saving the segmentation ## \n')
# seg_s3d_array = np.swapaxes(seg_s3d_array,0,2)
# # save_itk(seg_s3d_array, origin, spacings, path2points+'/seg_s3d.nrrd')
# save_itk_keeping_header(new_image=seg_s3d_array, original_image=seg_array_good_header, filename=path2points+'/seg_s3d.nrrd')

# print(" ## Pulmonary artery wall: Saved segmentation with pulmonary artery wall pushed inside ## \n")

# ----------------------------------------------------------------------------------------------
# Open the aorta and pulmonary artery
# ----------------------------------------------------------------------------------------------
# print(' ## Step 4/10: Opening arteries ## \n')
# seg_s3d_array, header = nrrd.read(path2points+'seg_s3d.nrrd')

# aorta_slicer_nrrd = path2points+'/aorta_slicer.nrrd'
# aorta_slicer_label = 0

# aorta_slicer_array, header = nrrd.read(aorta_slicer_nrrd)

# print(' ## Cropping major vessels: Slicing the aortic wall ## \n')
# seg_s3e_array = add_masks_replace_only(seg_s3d_array, aorta_slicer_array, aorta_slicer_label, Ao_wall_label)

# PArt_slicer_nrrd = path2points+'/PArt_slicer.nrrd'
# PArt_slicer_label = 0

# PArt_slicer_array, header = nrrd.read(PArt_slicer_nrrd)

# print(' ## Opening major vessels: Slicing the pulmonary artery wall ## \n')
# seg_s3e_array = add_masks_replace_only(seg_s3e_array, PArt_slicer_array, PArt_slicer_label, PArt_wall_label)

# cut_labels.open_artery(path2image      = f"{path2points}/seg_s3d.nrrd", 
# 					   myo_label       = Ao_wall_label, 
# 					   artery_bp_label = Ao_BP_label, 
# 					   ventricle_label = LV_BP_label, 
# 					   cut_ratio       = 0.95,
# 					   save_filename   = f"{path2points}/seg_s3d_aorta.nrrd")

# cut_labels.open_artery(path2image      = f"{path2points}/seg_s3d_aorta.nrrd", 
# 					   myo_label       = PArt_wall_label, 
# 					   artery_bp_label = PArt_BP_label, 
# 					   ventricle_label = RV_BP_label, 
# 					   cut_ratio       = 0.85,
# 					   save_filename   = f"{path2points}/seg_s3e.nrrd")

# ----------------------------------------------------------------------------------------------
# Format and save the segmentation
# ----------------------------------------------------------------------------------------------
# print(' ## Cropping major vessels: Formatting and saving the segmentation ## \n')
# seg_s3e_array = np.swapaxes(seg_s3e_array,0,2)
# save_itk(seg_s3e_array, origin, spacings, path2points+'/seg_s3e.nrrd')
# save_itk_keeping_header(new_image=seg_s3e_array, original_image=seg_array_good_header, filename=path2points+'/seg_s3e.nrrd')

# Sometimes the wall is over the slicer, so another connected appears with only myocardium, so we need to rerun the connected component.

## NOTE: run connected component keeping only the biggest.

# print(' ## Running connected component in the aorta wall## \n')
# seg_s3e_nrrd = f"{path2points}/seg_s3e.nrrd"
# seg_aorta_wall_cc_array = connected_component(seg_s3e_nrrd, Ao_tip_seed, Ao_wall_label,path2points)
# seg_aorta_wall_cc_array = np.swapaxes(seg_aorta_wall_cc_array,0,2)
# save_itk_keeping_header(new_image=seg_aorta_wall_cc_array, original_image=seg_array_good_header, filename=f"{path2points}/seg_aorta_wall_cc.nrrd")

# seg_aorta_wall_biggest_cc = connected_component_keep_biggest(imga_nrrd=f"{path2points}/tmp/CC.nrrd",
# 															 layer=Ao_wall_label)
# seg_aorta_wall_biggest_cc = np.swapaxes(seg_aorta_wall_biggest_cc,0,2)
# save_itk_keeping_header(new_image=seg_aorta_wall_biggest_cc, original_image=seg_array_good_header, filename=f"{path2points}/seg_aorta_wall_biggest_cc.nrrd")

# print(' ## Running connected component in the pulmonary artery wall ## \n')
# seg_aorta_wall_biggest_cc_nrrd = f"{path2points}/seg_aorta_wall_biggest_cc.nrrd"
# seg_PA_wall_cc_array = connected_component(seg_aorta_wall_biggest_cc_nrrd, PArt_tip_seed, PArt_wall_label,path2points)
# seg_PA_wall_cc_array = np.swapaxes(seg_PA_wall_cc_array,0,2)

# save_itk_keeping_header(new_image=seg_PA_wall_cc_array, original_image=seg_array_good_header, filename=f"{path2points}/seg_PA_wall_cc.nrrd")


# seg_PA_wall_biggest_cc = connected_component_keep_biggest(imga_nrrd=f"{path2points}/tmp/CC.nrrd",
# 															 layer=PArt_wall_label)
# seg_PA_wall_biggest_cc = np.swapaxes(seg_PA_wall_biggest_cc,0,2)
# save_itk_keeping_header(new_image=seg_PA_wall_biggest_cc, original_image=seg_array_good_header, filename=f"{path2points}/seg_PA_wall_biggest_cc.nrrd")

print(" ## Cropping major vessels: Saved segmentation with aorta and pulmonary artery sliced ## \n")


# # ----------------------------------------------------------------------------------------------
# # Prepare the seeds for the tips of the aorta and pulmonary artery
# # ----------------------------------------------------------------------------------------------
# Ao_tip_seed = points_data['Ao_tip']
# PArt_tip_seed = points_data['PArt_tip']
# Ao_wall_tip_seed = points_data['Ao_WT_tip']
# PArt_wall_tip_seed = points_data['PArt_WT_tip']

# # ----------------------------------------------------------------------------------------------
# # Crop the aorta and save the segmentation
# # ----------------------------------------------------------------------------------------------
# print(' ## Cropping major vessels: Cropping the aorta ## \n')
# seg_s3e_nrrd = path2points+'/seg_s3e.nrrd'
# seg_s3f_array = connected_component(seg_s3e_nrrd, Ao_tip_seed, Ao_BP_label,path2points)
# seg_s3f_array = np.swapaxes(seg_s3f_array,0,2)
# # save_itk(seg_s3f_array, origin, spacings, path2points+'/seg_s3f.nrrd')
# save_itk_keeping_header(new_image=seg_s3f_array, original_image=seg_array_good_header, filename=path2points+'/seg_s3f.nrrd')


# # ----------------------------------------------------------------------------------------------
# # Crop the pulmonary artery and save the segmentation
# # ----------------------------------------------------------------------------------------------
# print(' ## Cropping major vessels: Cropping the pulmonary artery ## \n')
# seg_s3f_nrrd = path2points+'/seg_s3f.nrrd'
# seg_s3f_array = connected_component(seg_s3f_nrrd, PArt_tip_seed, PArt_BP_label,path2points)
# seg_s3f_array = np.swapaxes(seg_s3f_array,0,2)
# # save_itk(seg_s3f_array, origin, spacings, path2points+'/seg_s3f.nrrd')
# save_itk_keeping_header(new_image=seg_s3f_array, original_image=seg_array_good_header, filename=path2points+'/seg_s3f.nrrd')

# ----------------------------------------------------------------------------------------------
# Create the RV myocardium
# ----------------------------------------------------------------------------------------------
print(' ## Step 5/10: Creating the right ventricular myocardium: ## \n')
print(' ## RV myo: Executing distance map ## \n')
RV_BP_DistMap = distance_map(path2points+'seg_s3c_PA.nrrd',RV_BP_label)
print(' ## RV myo: Writing temporary image ## \n')
sitk.WriteImage(RV_BP_DistMap,path2points+'/tmp/RV_BP_DistMap.nrrd',True)

print(' ## RV myo: Thresholding distance filter ## \n')
RV_myo = threshold_filter_nrrd(path2points+'/tmp/RV_BP_DistMap.nrrd',0,RV_WT)
sitk.WriteImage(RV_myo,path2points+'/tmp/RV_myo.nrrd',True)

print(' ## RV myo: Adding right ventricular myocardium to segmentation ## \n')
RV_myo_array, header = nrrd.read(path2points+'/tmp/RV_myo.nrrd')
seg_s3e_array, header = nrrd.read(path2points+'seg_s3c_PA.nrrd')
RV_myo_array = add_masks_replace(RV_myo_array,RV_myo_array,RV_myo_label)
seg_s3g_array = add_masks_replace_only(seg_s3e_array,RV_myo_array,RV_myo_label,Ao_wall_label) 

# ----------------------------------------------------------------------------------------------
# Format and save the segmentation
# ----------------------------------------------------------------------------------------------
print(' ## RV myo: Formatting and saving the segmentation ## \n')
seg_s3g_array = np.swapaxes(seg_s3g_array,0,2)
# save_itk(seg_s3g_array, origin, spacings, path2points+'/seg_s3g.nrrd')
save_itk_keeping_header(new_image=seg_s3g_array, original_image=seg_array_good_header, filename=path2points+'/seg_s3g.nrrd')

print(" ## RV myo: Saved segmentation with right ventricular myocardium added ## \n")

# ----------------------------------------------------------------------------------------------
# Create the LA myocardium
# ----------------------------------------------------------------------------------------------
print(' ## Step 6/10: Creating the left atrial myocardium: ## \n')
print(' ## LA myo: Executing distance map ## \n')
LA_BP_DistMap = distance_map(path2points+'seg_s3g.nrrd',LA_BP_label)
print(' ## LA myo: Writing temporary image ## \n')
sitk.WriteImage(LA_BP_DistMap,path2points+'/tmp/LA_BP_DistMap.nrrd',True)

print(' ## LA myo: Thresholding distance filter ## \n')
LA_myo = threshold_filter_nrrd(path2points+'/tmp/LA_BP_DistMap.nrrd',0,LA_WT)
sitk.WriteImage(LA_myo,path2points+'/tmp/LA_myo.nrrd',True)

print(' ## LA myo: left atrial myocardium to segmentation ## \n')
LA_myo_array, header = nrrd.read(path2points+'/tmp/LA_myo.nrrd')
seg_s3g_array, header = nrrd.read(path2points+'seg_s3g.nrrd')
LA_myo_array = add_masks_replace(LA_myo_array,LA_myo_array,LA_myo_label)
seg_s3h_array = add_masks_replace_only(seg_s3g_array,LA_myo_array,LA_myo_label,RA_BP_label)
seg_s3h_array = add_masks_replace_only(seg_s3h_array,LA_myo_array,LA_myo_label,SVC_label)

# ----------------------------------------------------------------------------------------------
# Format and save the segmentation
# ----------------------------------------------------------------------------------------------
print(' ## LA myo: Formatting and saving the segmentation ## \n')
seg_s3h_array = np.swapaxes(seg_s3h_array,0,2)
# save_itk(seg_s3h_array, origin, spacings, path2points+'/seg_s3h.nrrd')
save_itk_keeping_header(new_image=seg_s3h_array, original_image=seg_array_good_header, filename=path2points+'/seg_s3h.nrrd')

print(" ## LA myo: Saved segmentation with left atrial myocardium added ## \n")

# ----------------------------------------------------------------------------------------------
# Create the RA myocardium
# ----------------------------------------------------------------------------------------------
print(' ## Step 7/10: Creating the right atrial myocardium: ## \n')
print(' ## RA myo: Executing distance map ## \n')
RA_BP_DistMap = distance_map(path2points+'seg_s3h.nrrd',RA_BP_label)
print(' ## RA myo: Writing temporary image ## \n')
sitk.WriteImage(RA_BP_DistMap,path2points+'/tmp/RA_BP_DistMap.nrrd',True)

print(' ## RA myo: Thresholding distance filter ## \n')
RA_myo = threshold_filter_nrrd(path2points+'/tmp/RA_BP_DistMap.nrrd',0,RA_WT)
sitk.WriteImage(RA_myo,path2points+'/tmp/RA_myo.nrrd',True)
"""
print(' ## RA myo: Adding right atrial myocardium to segmentation ## \n')
RA_myo_array, header = nrrd.read(path2points+'/tmp/RA_myo.nrrd')
seg_s3hi_array, header = nrrd.read(path2points+'seg_s3h.nrrd')

RA_myo_array = add_masks_replace(RA_myo_array,RA_myo_array,RA_myo_label)
"""
seg_s3i_array = add_masks_replace_only(seg_s3hi_array,RA_myo_array,RA_myo_label,RPV1_label)

# ----------------------------------------------------------------------------------------------
# Format and save the segmentation
# ----------------------------------------------------------------------------------------------
print(' ## RA myo: Formatting and saving the segmentation ## \n')
seg_s3i_array = np.swapaxes(seg_s3i_array,0,2)
# save_itk(seg_s3i_array, origin, spacings, path2points+'/seg_s3i.nrrd')
save_itk_keeping_header(new_image=seg_s3i_array, original_image=seg_array_good_header, filename=path2points+'/seg_s3i.nrrd')

print(" ## RA myo: Saved segmentation with right atrial myocardium added ## \n")

# ----------------------------------------------------------------------------------------------
# Right atrium is pushed by the left atrium
# ----------------------------------------------------------------------------------------------
print(' ## RA myo: Pushing the right atrium with the left atrium ## \n')
seg_s3j_array = push_inside(path2points,path2points+'seg_s3i.nrrd',LA_myo_label,RA_myo_label,RA_BP_label,RA_WT)

# ----------------------------------------------------------------------------------------------
# Format and save the segmentation
# ----------------------------------------------------------------------------------------------
print(' ## RA myo: Formatting and saving the segmentation ## \n')
seg_s3j_array = np.swapaxes(seg_s3j_array,0,2)
# save_itk(seg_s3j_array, origin, spacings, path2points+'/seg_s3j.nrrd')
save_itk_keeping_header(new_image=seg_s3j_array, original_image=seg_array_good_header, filename=path2points+'/seg_s3j.nrrd')

print(" ## RA myo: Saved segmentation with right atrium pushed by the left atrium ## \n")

# ----------------------------------------------------------------------------------------------
# Right atrium is pushed by the aorta
# ----------------------------------------------------------------------------------------------
print(' ## RA myo: Pushing the right atrium with the aorta ## \n')
seg_s3k_array = push_inside(path2points,path2points+'seg_s3j.nrrd',Ao_wall_label,RA_myo_label,RA_BP_label,RA_WT)

print(' ## RA myo: Formatting and saving the segmentation ## \n')
seg_s3k_array = np.swapaxes(seg_s3k_array,0,2)
# save_itk(seg_s3k_array, origin, spacings, path2points+'/seg_s3k.nrrd')
save_itk_keeping_header(new_image=seg_s3k_array, original_image=seg_array_good_header, filename=path2points+'/seg_s3k.nrrd')

print(" ## RA myo: Saved segmentation with right atrium pushed by the aorta ## \n")



print(' ## SVC: Pushing the SVC with the aorta to have space for the SVC ring## \n')

# Step 1: Create disk in SVC hole with RA myo thickness

# We had the RA_myo_array from when we created the RA myo
SVC_disk = and_filter(imga_array=seg_s3k_array, imgb_array=RA_myo_array, label_a=SVC_label, new_label=SVC_label+100)

# Step 2: Intersect disk with the dilated aorta using the ring thickness and assign it to RA myo

dilated_aorta_image = threshold_filter_nrrd(path2points+'/tmp/Ao_DistMap.nrrd',0,Ao_WT + rings_thickness) # The distance map is from the blood pool
sitk.WriteImage(dilated_aorta_image,path2points+'/tmp/aorta_distmap_svc_ring.nrrd',True)
dilated_aorta, _ = nrrd.read(path2points+'/tmp/aorta_distmap_svc_ring.nrrd')  # Now it's an array


intersected_disk_array = and_filter(SVC_disk, dilated_aorta, SVC_label+100, RA_myo_label) # This is a semilunar disk of RA myo

# Same steps as always to replace a mask
intersected_disk_array = add_masks_replace(intersected_disk_array,intersected_disk_array,RA_myo_label)
seg_s3k_disk_array = add_masks_replace_only(seg_s3k_array,intersected_disk_array,RA_myo_label,SVC_label)

# Step 3: Intersect SVC with the dilated aorta using the ring thickness and set to 0

intersected_SVC_array = and_filter(seg_s3k_disk_array, dilated_aorta, SVC_label, 0) # We ate away a piece of the SVC cylinder with the aorta

# Same steps as always to replace a mask
intersected_SVC_array = add_masks_replace(intersected_SVC_array,intersected_SVC_array,0)
seg_s3k_array = add_masks_replace_only(seg_s3k_disk_array,intersected_SVC_array,0,SVC_label)

# seg_s3k_array = push_inside(path2points,path2points+'seg_s3k.nrrd',Ao_wall_label,RA_myo_label,SVC_label,RA_WT)

# ----------------------------------------------------------------------------------------------
# Format and save the segmentation
# ----------------------------------------------------------------------------------------------
print(' ## RA myo: Formatting and saving the segmentation ## \n')
# DO NOT SWAP AXES THIS TIME!!! Otherwise the heart will be mirrored. An alternative solution is to load the array from file after saving it, instead of reusing the array.

# seg_s3k_array = np.swapaxes(seg_s3k_array,0,2)
# save_itk(seg_s3k_array, origin, spacings, path2points+'/seg_s3k.nrrd')
save_itk_keeping_header(new_image=seg_s3k_array, original_image=seg_array_good_header, filename=path2points+'/seg_s3k.nrrd')

print(" ## RA myo: Saved segmentation with right atrium pushed by the aorta ## \n")

# ----------------------------------------------------------------------------------------------
# Right atrium is pushed by the left ventricle
# ----------------------------------------------------------------------------------------------
print(' ## RA myo: Pushing the right atrium with the left ventricle ## \n')
seg_s3l_array = push_inside(path2points,path2points+'seg_s3k.nrrd',LV_myo_label,RA_myo_label,RA_BP_label,RA_WT)

# ----------------------------------------------------------------------------------------------
# Format and save the segmentation
# ----------------------------------------------------------------------------------------------
print(' ## RA myo: Formatting and saving the segmentation ## \n')
seg_s3l_array = np.swapaxes(seg_s3l_array,0,2)
# save_itk(seg_s3l_array, origin, spacings, path2points+'/seg_s3l.nrrd')
save_itk_keeping_header(new_image=seg_s3l_array, original_image=seg_array_good_header, filename=path2points+'/seg_s3l.nrrd')

print(" ## RA myo: Saved segmentation with right atrium pushed by the left ventricle ## \n")

# ----------------------------------------------------------------------------------------------
# Left atrium is pushed by the aorta
# ----------------------------------------------------------------------------------------------
print(' ## Step 8/10: LA myo: Pushing the left atrium with the aorta ## \n')
seg_s3m_array = push_inside(path2points,path2points+'seg_s3l.nrrd',Ao_wall_label,LA_myo_label,LA_BP_label,LA_WT)

# ----------------------------------------------------------------------------------------------
# Format and save the segmentation
# ----------------------------------------------------------------------------------------------
print(' ## LA myo: Formatting and saving the segmentation ## \n')
seg_s3m_array = np.swapaxes(seg_s3m_array,0,2)
# save_itk(seg_s3m_array, origin, spacings, path2points+'/seg_s3m.nrrd')
save_itk_keeping_header(new_image=seg_s3m_array, original_image=seg_array_good_header, filename=path2points+'/seg_s3m.nrrd')

print(" ## LA myo: Saved segmentation with left atrium pushed by the aorta ## \n")

# ----------------------------------------------------------------------------------------------
# Pulmonary artery is pushed by the aorta
# ----------------------------------------------------------------------------------------------
print(' ## Step 9/10: PArt wall: Pushing the pulmonary artery with the aorta ## \n')
seg_s3n_array = push_inside(path2points,path2points+'seg_s3m.nrrd',Ao_wall_label,PArt_wall_label,PArt_BP_label,PArt_WT)

# ----------------------------------------------------------------------------------------------
# Format and save the segmentation
# ----------------------------------------------------------------------------------------------
print(' ## PArt wall: Formatting and saving the segmentation ## \n')
seg_s3n_array = np.swapaxes(seg_s3n_array,0,2)
# save_itk(seg_s3n_array, origin, spacings, path2points+'/seg_s3n.nrrd')
save_itk_keeping_header(new_image=seg_s3n_array, original_image=seg_array_good_header, filename=path2points+'/seg_s3n.nrrd')

print(" ## PArt wall: Saved segmentation with pulmonary artery pushed by the aorta ## \n")

# ----------------------------------------------------------------------------------------------
# Pulmonary artery is pushed by the left ventricle
# ----------------------------------------------------------------------------------------------
print(' ## PArt wall: Pushing the pulmonary artery with the left ventricle ## \n')
seg_s3o_array = push_inside(path2points,path2points+'seg_s3n.nrrd',LV_myo_label,PArt_wall_label,PArt_BP_label,PArt_WT)

# ----------------------------------------------------------------------------------------------
# Format and save the segmentation
# ----------------------------------------------------------------------------------------------
print(' ## PArt wall: Formatting and saving the segmentation ## \n')
seg_s3o_array = np.swapaxes(seg_s3o_array,0,2)
# save_itk(seg_s3o_array, origin, spacings, path2points+'/seg_s3o.nrrd')
save_itk_keeping_header(new_image=seg_s3o_array, original_image=seg_array_good_header, filename=path2points+'/seg_s3o.nrrd')

print(" ## PArt wall: Saved segmentation pulmonary artery pushed by the left ventricle ## \n")

# ----------------------------------------------------------------------------------------------
# Right ventricle is pushed by the aorta
# ----------------------------------------------------------------------------------------------
print(' ## Step 10/10: RV myo: Pushing the right ventricle with the aorta ## \n')
seg_s3p_array = push_inside(path2points,path2points+'seg_s3o.nrrd',Ao_wall_label,RV_myo_label,RV_BP_label,RV_WT)

# ----------------------------------------------------------------------------------------------
# Format and save the segmentation
# ----------------------------------------------------------------------------------------------
print(' ## RV myo: Formatting and saving the segmentation ## \n')
seg_s3p_array = np.swapaxes(seg_s3p_array,0,2)
# save_itk(seg_s3p_array, origin, spacings, path2points+'/seg_s3p.nrrd')
save_itk_keeping_header(new_image=seg_s3p_array, original_image=seg_array_good_header, filename=path2points+'/seg_s3p.nrrd')

print(" ## RV myo: Saved segmentation with right ventricle pushed by the aorta ## \n")

print("Things to check: The aorta and the PA have only one piece each.")