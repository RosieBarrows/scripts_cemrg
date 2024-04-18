# 4ch pipeline - PART I

*This first half of the pipeline takes the segmentation as an input and performs a series of post-processing steps, including creating the myocardium, valve planes and vein rings*

Instructions for using four-chamber mesh generation pipeline

**Set-up**

To set up conda environment to run scripts:

conda create -n py4ch python=3.9

pip install:    
        numpy~=1.23.0 \
        meshio~=5.3.4 \
        vtk~=9.2.6 \
        tqdm~=4.65.0 \
        pyvista~=0.39.1 \
        trimesh~=3.22.1 \
        scipy~=1.10.1  \

conda deactivate
conda activate py4ch

**Segmentation**

1) Make folder (e.g. heart_folder_) then make folder “segmentations” inside
2) From /scripts_cemrgapp/h_template/segmentations, copy “points.txt”, “labels.txt”, “origin_spacing.txt” and “origin_spacing_labels.txt” templates into "segmentations" folder
3) Move end-diastolic image CT image into heart_folder/segmentations folder and call that directory "ct"
4) Move segmented CT into "segmentations" folder.
5) In Seg3D (or alternative), manually split and relabel pulmonary veins:
6) Export the segmentation as “seg_corrected.nrrd” with the following tags:
        LV blood pool = 1
        LV myocardium = 2
        RV blood pool = 3
        LA blood pool = 4
        RA blood pool = 5
        aorta blood pool = 6
        pulmonary artery blood pool = 7
        left superior pulmonary vein = 8
        left inferior pulmonary vein = 9
        right superior pulmonary vein = 10
        right inferior pulmonary vein = 11
        left atrial appendage = 12

(hard-coded tags replaced by .json files in CemrgApp implementation)

7) From seg_scripts,
   ```
   python find_origin_and_spacing.py [/heart_folder/segmentations]
   ```
9) Copy origin and spacing into “origin_spacing.txt” file - MAKE SURE THERE ARE NO COMMAS/QUOTE MARKS/ETC

10) Select 3 points for SVC and IVC cylinders (pixel) and save in “points.txt”
11) Select 3 points for Ao_slicer and for PArt_slicer and save in “points.txt”
    a. Note: be careful not to crop the PArt too close to the RV_BP
12) From seg_scripts,
    ```
    python create_cylinders.py [/heart_folder/segmentations]
    ```
13) From seg_scripts,
    ```
    python create_svc_ivc.py [/heart_folder/segmentations]
    ```

14) Select:
    3 points for SVC_slicer 
    1 point for SVC_tip
    3 points for IVC_slicer
    1 point for IVC_tip
    1 point for Ao_tip
    1 point for PArt_tip

    and save in “points.txt”
    
    a. Note: SVC/IVC tips are used for removing sections protruding from wrong side of RA

15) From seg_scripts,
    ```
    python create_slicers.py [/heart_folder/segmentations]
    ```

17) From seg_scripts,
    ```
    python crop_svc_ivc.py [/heart_folder/segmentations]
    ```

18) From seg_scripts.
    ```
    python create_myo.py [/heart_folder/segmentations]
    ```

19) Select a seed point for the Ao_wall_tip and PArt_wall_tip
20) From seg_scripts,
    ```
    python create_valve_planes.py [/heart_folder/segmentations]
    ```

21) Load the most recent .nrrd and inspect the segmentation. From the clean_seg.py scripts, leave uncommented those modifications that would improve your segmentation (it is best that the rest stay commented out).

    From seg_scripts,
    ```
    python clean_seg.py [/heart_folder/segmentations]
    ```

20) Load the most recent .nrrd file and inspect the segmentation. Any final manual corrects can be completed at this point. 
21) Export segmentation as “seg_final.nrrd”

22) From seg_scripts,
    ```
    bash segSmoothing.sh [/heart_folder/segmentations]
    ```

**Meshing**

1) From meshing,
   ```
   meshing.py [/heart_folder]
   ```

2) ```cd [/heart_folder/meshing/]```
3) ```bash generate_heart_mesh_Cesare.sh```

5) From meshing,
   ```
   bash mesh_post.sh [/heart_folder]
   ```
   
   **watch out mesh tool path currently hard-coded **














