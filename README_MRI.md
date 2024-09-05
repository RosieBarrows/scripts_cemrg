Instructions for using four-chamber mesh generation pipeline
=========================================================================================================

Segmentation
=========================================================================================================
1) Make folder (e.g. henry) in /data/Dropbox then make folder `segmentations' inside
<!-- 2) Copy `points.txt`, `labels.txt`, `origin_spacing.txt` and `origin_spacing_labels.txt` templates into `segmentations` folder (currently in `h_templates` folder) -->
5) Manually split and relabel pulmonary veins:
6) Export the segmentation as `seg_corrected.nrrd` with the new tags:
- LSPV = 8
- LIPV = 9
- RSPV = 10
- RIPV = 11
- LAA = 12
1) Make sure the image is in a suitable orientation (by opening it in ITKsnap, for instance). If not, you can use `c3d seg_corrected_badly_oriented.nii -orient RAI -o seg_corrected_RAI.nii.gz` to convert it to RAI orientation.
5) Resample to 0.3 mm isotropic using the `img.py` script. The final segmentation name has to be `seg_corrected.nrrd` for the following steps to work.

6) Run `create_origin_spacing_files` in `img.py`. Otherwise, select manually in itksnap the origin and spacing and add it to origin_spacing.txt

	To select points in itksnap:
	
	1) Load the segmentation as the main image. If it doesn't work, make sure that your image has the right orientation (step 1). 
	In `Tools` > `Image information...` you'll be able to see the information about the spacing and the origin.
	2) To select points in itksnap, once the main image has been loaded, drag the segmentation again, this time `Load as segmentation`. 
	3) Click `update` in the 3D window to see the 3D representation.
	4) On the left bar, at the bottom, you will see a `3D Toolbar`. Click the second icon (looks like a cross-hair with a cursor on top).
	5) Click on any window (including the 3D view) to select a point. You can refine in the 2D windows by placing the cursor on one of them and scrolling. Your selected point will be (in voxel index) at the left bar, under `Cursor position (x,y,z)`. That is the information you will need to copy.
	6) To manipulate the geometry without selecting points you need to click the first icon under the `3D Toolbar` (looks like a 3D axis with a curved arrow).

	This process is almost automated in CEMRGApp but it's not stable. Ask Jose about it if interested.

	If the scripts complain about a missing .txt file check the templates in `/scripts_cemrg/h_template/segmentations`.

9) Select 3 points for SVC and IVC cylinders (pixel) and save in `points.txt`

11) Run `create_cylinders.py /heart_folder/segmentations`. Note that if you visualize the SVC and IVC at this stage it will look off. Just continue with the pipeline and check them in the following step.

12) Run `create_svc_ivc.py /heart_folder/segmentations`. At this point the segmentation should have the venae cavae incorporated.

1) Copy the file `vc_joint.json` from `h_template/segmentations` to your segementation folder. The information there will be used to add a piece of the RA to the venae cavae to avoid areas of too few cotact. Ideally, you just want to add enough so there are no sharp corners. The values go from 0 to 1, where 1 means that 100% of the cylinder for the vena cava is added to RA. Run `cut_vessels.py /heart_folder/segmentations`. Check that the venae cava have some base, otherwise adjust the threhsold values. 

17) Run `create_myo.py /heart_folder/segmentations`

19) Run `create_valve_planes.py /heart_folder/segmentations`

20) Run `clean_seg.py /heart_folder/segmentations`

21) Run `correct_pulmonary_valve.py /heart_folder/segmentations`

21) Check for any lose voxels. Make sure all the chambers are closed by the valve planes.
22) Export segmentation as `seg_final.nrrd`

At this point, talk with Cristobal for the next steps.

<!-- 
23) Run `prepare_for_segsmooth.py`

23) Run `segSmoothing.sh /heart_folder/segmentations`. If you get a warning about "ghost labels", run it as `segSmoothing.sh /heart_folder/segmentations 2`. 

24) Run `align_segmentation.py /heart_folder/segmentations`. Not very well tested, check the output of this compared to the dicom/nifti images or with the first segmentation of the pipeline.

23) Check all the labels one last time. It is common that lose voxels appear, so run connected components in all the labels. Export it as `seg_final_smooth_corrected.nrrd`. Convert it to `.inr` using `segconvert --resolution 0.15 0.15 0.15` -->

Meshing
=========================================================================================================
1) Run `python3 meshing.py /heart_folder`

2) cd into the meshing folder (inside heart folder)
3) Run `bash generate_heart_mesh_Cesare.sh`

4) Run `bash mesh_post.sh /heart_folder` ## watch out mesh tool path currently hard-coded ##














