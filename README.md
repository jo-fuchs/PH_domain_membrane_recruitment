# Membrane recruitment analysis

![](/MemRecruitment_Method.jpg?raw=true "Example of membrane analysis workflow")

## Membrane recruitment analysis part 1
ImageJ-Macro: Save_profiles.ijm

Run for every image separately after:
 - opening Image-Stack with 4 channels
 - Lines stored in ROI-manager

Goals:
 - check &/or create LineProfiles-folder
 - Create line profile of a stored line in the ROI manager in all 4 channels
 - Scale Y-axis to µm instead of pixeles
 - Save final table for each ROI with the name: "ImageName_RoiName"
 - Save ROIs


## Membrane recruitment analysis part 2
R-script Membrane_recruitment.R with function merge_lines()

Goals:
 - Normalizing line profiles to peaks of cell membrane marker (here: F-actin) to create merge
 - Read out Membrane to Cytosol ratio of Intensities
 - Create quality control graphs and summary Figures

