Instructions for using four-chamber mesh generation pipeline
=========================================================================================================

Segmentation
=========================================================================================================
1) Make folder (e.g. henry) in /data/Dropbox then make folder “segmentations” inside
2) Copy “points.txt”, “labels.txt”, “origin_spacing.txt” and “origin_spacing_labels.txt” templates into "segmentations" folder (currently in “h_templates” folder)
5) Manually split and relabel pulmonary veins:
6) Export the segmentation as “seg_corrected.nrrd” with the new tags:
        LSPV = 8
        LIPV = 9
        RSPV = 10
        RIPV = 11
        LAA = 12
 Resample to 0.3 mm isotropic using the `img.py` script.

7) Run select manually in itksnap the origin and spacing and add it to origin_spacing.txt

9) Select 3 points for SVC and IVC cylinders (pixel) and save in “points.txt”
10) Select 3 points for Ao_slicer and for PArt_slicer and save in “points.txt”
    a. Note: be careful not to crop the PArt too close to the RV_BP
11) Run “create_cylinders.py /heart_folder/segmentations”

12) Run “create_svc_ivc.py /heart_folder/segmentations”. At this point the segmentation should have the venae cavae incorporated.

13) Select:
    3 points for SVC_slicer 
    1 point for SVC_tip
    3 points for IVC_slicer
    1 point for IVC_tip
    1 point for Ao_tip
    1 point for PArt_tip

    and save in “points.txt”
    
    a. Note: SVC/IVC tips are used for removing sections protruding from wrong side of RA

14) Run "create_slicers.py /heart_folder/segmentations/" 
16) Run “crop_svc_ivc.py /heart_folder/segmentations”

If you need to re-do the slicers of the aorta or pulmonary artery, you have to re run 1_create_cylinders.

17) Run “create_myo.py /heart_folder/segmentations”

19) Run “create_valve_planes.py /heart_folder/segmentations”

20) Run “clean_seg.py /heart_folder/segmentations”

21) Load .nrrd and manually clean valve ring
22) Export segmentation as “seg_final.nrrd”

23) Run “segSmoothing.sh /heart_folder/segmentations”


Meshing
=========================================================================================================
1) Run “meshing.py /heart_folder”

2) cd into the meshing folder (inside heart folder)
3) Run “bash generate_heart_mesh_Cesare.sh”

4) Run “bash mesh_post.sh /heart_folder” ## watch out mesh tool path currently hard-coded ##













