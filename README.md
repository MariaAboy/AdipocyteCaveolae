# Adipocyte caveolae spatial description
R scripts and ImageJ macros developed for the semi-automated detection of caveolae in EM microscopy images of adipocytes

## CompleteCaveolaeAnalysis.ijm
This ImgeJ macro takes previously detected caveolae coordinates using the CellCounter plugin and manually delineated PM and LD, and calculates:
1) Through caveolae PM/LD distance
2) Global PM/LD distance
3) Cytoplasmic coordinates

## TEMscript.R
This script uses the caveolar and cytoplasmic coordinates. It tessellates the cytoplasmic area applying a square grid of 20x20 pixels (29.6x29.6 nm). Manually detected caveolae coordinates are projected to the closer grid points. Global cellular PM-LD distances are calculated as the minimum euclidian distance of each grid point to the PM and the LD. Through-caveolae PM-LD distance was calculated as the PM-LD distance of the grid point closer to the manually assigned caveolae coordinate.

## ITEMscript.R
This script works similarly to TEMscript.R but it is adapted to inlcude gold label coordiantes in Immunogold images

## TIRFvesicles.ijm
This ImageJ macro was designed to work with image stacks of two fluorescence channels + brightfield, and an additional channel of TIRM. It will measure cell morphology, epifluorescene intensity of BODIPY, and Cav1, and TIRFM vesicle number and intensity. It follows the following steps:
1) creates a stack with the four channels
2) Asks user to draw cells and add to ROI manager in Cav1 channel
3) Measures cell morphology, Cav1 intensity, Bodipy intensity, detects vesicles, counts vesicles and measures vesicle intensity for each cell
4) Asks user to draw main LD
5) Measures LD size

The Macro outputs 12 files:
1) Stack image of 4 channels and mask of detected vesicles
2) 3 zip files for ROIs of cells, vesicles and LDs
3) files with suffix BODIPY, CAVEOLIN, VESICLES contain cell morphology measures and intensity measures for the BODIPY, Cav1 and TIRF channel for each cell
4) file with suffix VESICLECOUNT contains counts of vesicles for each cell
5) file with suffix LIPIDDROPLETS contains area of the meain LD of each cell
6) file with suffix SUMMARY contains vesicle counts, cell morphology measures, LD area and intensity measures for BODIPY, Cav1 and TIRF channel for all cells in the same table
7) file with suffix VESICLEMEASURE contains vesicle coordinates and intensity measures for each vesicle of each cell

