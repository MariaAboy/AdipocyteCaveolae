//This macro will measure cell morphology, Bodipy, Cav1 and TIRF vesicle number and intensity
//The procedure is SEMI-automated, and will need user input
//The input consists of two images:
//One image is a stack of three channels: Bodipy, Cav1 and brightfield
//The other image is a TIRF image of Cav1 vesicles
//What the macro does:
//1-creates a stack with the four channels
//2-Asks user to draw cells and add to ROI manager in Cav1 channel
//3-Measures cell morphology, Cav1 intensity, Bodipy intensity, detects vesicles, counts vesicles and measures vesicle intensity for each cell
//4-Asks user to draw main LD
//5-Measures LD size
//Macro outputs 12 files:
//Stack image of 4 channels and mask of detected vesicles
//3 zip files for ROIs of cells, vesicles and LDs
//files with suffix BODIPY, CAVEOLIN, VESICLES contain cell morphology measures and intensity measures for the Bodipy, Cav1 and TIRF channel for each cell
//file with suffix VESICLECOUNT contains counts of vesicles for each cell
//file with suffix LIPIDDROPLETS contains area of the meain LD of each cell
//file with suffix SUMMARY contains vesicle counts, cell morphology measures, LD area and intensity measures for bodipy, Cav1 and TIRF channel for all cells in the same table
//file with suffix VESICLEMEASURE contains vesicle coordinates and intensity measures for each vesicle of each cell

//macro writes output files in same folder of input files
//when saving, in case of matching filenames, the macro will overwrite without notice
//Some important comments within code below, please read


selectImage(1);
title1 = getTitle;
dir = getDirectory("image"); //both input images in same directory
getDimensions(width1, height1, channels1, slices1, frames1);
selectImage(2);
title2 = getTitle;
getDimensions(width2, height2, channels2, slices2, frames2);

selectImage(1);
if (channels1 == 3 ) {
	run("Add Slice", "add=channel");
	selectImage(title2);
	run("Copy");
	selectWindow(title1);
	Stack.setChannel(2) ;
	run("Paste");
	run("Red");
	selectWindow(title2);
	run("Close");
	selectWindow(title1);
	saveAs("Tiff", dir+ "ZZ" + title1 + "_STACK" + ".tiff");
			
			
	} else { selectImage(2);
		run("Add Slice", "add=channel");
		selectImage(title1);
		run("Copy");
		selectWindow(title2);
		Stack.setChannel(2) ;
		run("Paste");
		run("Red");
		selectWindow(title1);
		run("Close");
		selectWindow(title2);
		saveAs("Tiff",  dir+"ZZ" + title2 + "_STACK" + ".tiff");
		
		}


title = getTitle;
dir = getDirectory("image");
run("Set Measurements...", "area mean min perimeter shape integrated redirect=None decimal=9");
selectImage(title);
Stack.setChannel(3) ;
waitForUser("Draw cells and add to ROI manager"); //you can press "T" to add a cell to the ROI manager after drawing the cell. Use polygon tool to draw cell
roiManager("deselect");
ncells=roiManager("count");
Stack.setChannel(3) ;
for (i=0; i<ncells; i++) {
	roiManager("select", i);
	Stack.setChannel(3) ;
	run("Measure");
}
							/////ARRAYS
							area=newArray(ncells);
							for (i=0; i<nResults; i++){
							  area[i] = getResult("Area",i);
							}
							mediaCav=newArray(ncells);
							for (i=0; i<nResults; i++){
							  mediaCav[i] = getResult("Mean",i);
							}
							minCav=newArray(ncells);
							for (i=0; i<nResults; i++){
							  minCav[i] = getResult("Min",i);
							}
							maxCav=newArray(ncells);
							for (i=0; i<nResults; i++){
							  maxCav[i] = getResult("Max",i);
							}
							per=newArray(ncells);
							for (i=0; i<nResults; i++){
							  per[i] = getResult("Perim.",i);
							}
							circ=newArray(ncells);
							for (i=0; i<nResults; i++){
							  circ[i] = getResult("Circ.",i);
							}
							InDenCav=newArray(ncells);
							for (i=0; i<nResults; i++){
							  InDenCav[i] = getResult("IntDen",i);
							}
							RInDenCav=newArray(ncells);
							for (i=0; i<nResults; i++){
							  RInDenCav[i] = getResult("RawIntDen",i);
							}
							AR=newArray(ncells);
							for (i=0; i<nResults; i++){
							  AR[i] = getResult("AR",i);
							}
							Roun=newArray(ncells);
							for (i=0; i<nResults; i++){
							  Roun[i] = getResult("Round",i);
							}
							sol=newArray(ncells);
							for (i=0; i<nResults; i++){
							 sol[i] = getResult("Solidity",i);
							}
												
							/////ARRAYS


selectWindow("Results");
saveAs("txt", dir + "zz_" + title + "_CAVEOLIN.txt");
selectWindow("Results");
run("Close");


selectImage(title);
Stack.setChannel(2) ;
roiManager("deselect");
for (i=0; i<ncells; i++) {
	roiManager("select", i);
	Stack.setChannel(2) ;
	run("Measure");
}
							/////ARRAYS
							mediaVes=newArray(ncells);
							for (i=0; i<nResults; i++){
							  mediaVes[i] = getResult("Mean",i);
							}
							minVes=newArray(ncells);
							for (i=0; i<nResults; i++){
							  minVes[i] = getResult("Min",i);
							}
							maxVes=newArray(ncells);
							for (i=0; i<nResults; i++){
							  maxVes[i] = getResult("Max",i);
							}
							InDenVes=newArray(ncells);
							for (i=0; i<nResults; i++){
							  InDenVes[i] = getResult("IntDen",i);
							}
							RInDenVes=newArray(ncells);
							for (i=0; i<nResults; i++){
							  RInDenVes[i] = getResult("RawIntDen",i);
							}
																			
							/////ARRAYS
selectWindow("Results");
saveAs("txt", dir + "zz_" + title + "_VESICLES.txt");
selectWindow("Results");
run("Close");

selectImage(title);
Stack.setChannel(1) ;
roiManager("deselect");
for (i=0; i<ncells; i++) {
	roiManager("select", i);
	Stack.setChannel(1) ;
	run("Measure");
}
							/////ARRAYS
							mediaBdp=newArray(ncells);
							for (i=0; i<nResults; i++){
							  mediaBdp[i] = getResult("Mean",i);
							}
							minBdp=newArray(ncells);
							for (i=0; i<nResults; i++){
							  minBdp[i] = getResult("Min",i);
							}
							maxBdp=newArray(ncells);
							for (i=0; i<nResults; i++){
							  maxBdp[i] = getResult("Max",i);
							}
							InDenBdp=newArray(ncells);
							for (i=0; i<nResults; i++){
							  InDenBdp[i] = getResult("IntDen",i);
							}
							RInDenBdp=newArray(ncells);
							for (i=0; i<nResults; i++){
							  RInDenBdp[i] = getResult("RawIntDen",i);
							}
																			
							/////ARRAYS
selectWindow("Results");
saveAs("txt", dir + "zz_" + title + "_BODIPY.txt");
selectWindow("Results");
run("Close");

roiManager("deselect");
run("Select None");

selectImage(title);
Stack.setChannel(2) ;
run("Find Maxima...", "noise=140 output=[Point Selection] exclude");
//roiManager("Add");


run("Create Mask");
//run("Invert"); //You may need to run this line if your mask is white

Mask = getTitle();




ncells=roiManager("count");
for (i=0; i<ncells; i++) {
	selectWindow(Mask);
	run("Duplicate...", " ");
	MaskDup=getTitle();
	roiManager("select", i);
	wait(200);  //we use this wait to clearly see what is happening. Delete this line for greater speed
	run("Clear Outside");
	wait(200); //we use this wait to clearly see what is happening. Delete this line for greater speed
	run("Select None");
	selectWindow(MaskDup);
	run("Find Maxima...", "noise=0 output=Count");
	wait(200);
	selectWindow(MaskDup);
	run("Close");
}	
							/////ARRAYS
							CountVes=newArray(ncells);
							for (i=0; i<nResults; i++){
							  CountVes[i] = getResult("Count",i);
							}															
							/////ARRAYS
							
roiManager("save", dir + "zz_" + title + "_CELLS.zip");
selectWindow("Results");
	
wait(200);
saveAs("txt", dir + "zz_" + title + "_VESICLECOUNT.txt");
selectWindow("Mask");
saveAs("Tiff", dir+ "zz_" + title + "_MASK" + ".tiff");
selectWindow("Results");
run("Close");


	//new
	for (i=0; i<ncells; i++) {
	run("Set Measurements...", "area mean min display redirect=None decimal=9");
	selectWindow("zz_" + title + "_MASK" + ".tiff");
	run("Duplicate...", " ");
	MaskDup=getTitle();
	roiManager("select", i);
	wait(200); //we use this wait to clearly see what is happening. Delete this line for greater speed
	run("Clear Outside");
	wait(200); //we use this wait to clearly see what is happening. Delete this line for greater speed
	run("Select None");
	selectWindow(MaskDup);
	run("Find Maxima...", "noise=0 output=[Point Selection]");
	roiManager("Add");
	selectWindow(title);
	Stack.setChannel(2) ;
	roiManager("select", ncells + i );
	roiManager("Measure");
	selectWindow(MaskDup);
	run("Close");
}
	

roiManager("save", dir + "zz_" + title + "_CELLS+VESICLES.zip");
selectWindow("Results");
saveAs("txt", dir + "zz_" + title + "_VESICLEMEASURE.txt");
wait(200);
selectWindow("Results");
run("Close");



selectWindow("zz_" + title + "_MASK" + ".tiff");
run("Close");

selectImage(title);

run("Set Measurements...", "area redirect=None decimal=9");
n2cells=roiManager("count");
for (i=0; i<ncells; i++) {
	roiManager("select", i);
	Stack.setChannel(1);
	waitForUser("Draw main LD of this cell");
	roiManager("add");
	roiManager("select", n2cells+i);
	roiManager("measure");
	
}	
							/////ARRAYS
							areaLD=newArray(ncells);
							for (i=0; i<nResults; i++){
							  areaLD[i] = getResult("Area",i);
							}															
							/////ARRAYS

selectWindow("Results");
saveAs("txt", dir + "zz_" + title + "_LIPIDDROPLETS.txt");
selectWindow("Results");
run("Close");

roiManager("save", dir + "zz_" + title + "_CELLSandLIPIDDROPLETS.zip");
selectWindow("ROI Manager");
run("Close");


	Array.show(area, per, circ, AR, Roun, sol, CountVes, areaLD, mediaCav, minCav, maxCav, InDenCav, RInDenCav, mediaVes, minVes, maxVes, InDenVes, RInDenVes, mediaBdp, minBdp, maxBdp, InDenBdp, RInDenBdp  );
	selectWindow("Arrays");
	saveAs("Results", dir + "zz_" + title +  "_SUMMARY" + ".txt");

	
run("Close All");
	selectWindow("zz_" + title +  "_SUMMARY" + ".txt");
	run("Close");