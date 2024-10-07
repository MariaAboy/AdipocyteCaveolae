//This macro analyses PM-LD distances through caveolae and globally
//It requires that the stitched EM image is open,
//caveolae manually counted previously using Cell Counter
//and ROIs for PM and LD previously delimited

//macro is semi-automated. It has several wait-for-user steps that require the user to follow instructions to load
//previously generated caveolae coordinates and PM and LD to ROI Manager


/////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
///////PART1 GET CAVEOLAE COORDINATES//////////////
////////////////////////////////////////////////////////


run("Set Scale...", "distance=676.4536 known=1 pixel=1 unit=µm global");

title = getTitle;
dir = getDirectory("image");

Dialog.create("Imagen");
Dialog.addString("Célula","");
Dialog.show();
cell=Dialog.getString;
Dialog.create("Imagen");
Dialog.addString("Parte","");
Dialog.show();
part=Dialog.getString;

//we get caveolae info
//detected caveolae coordinates will be passed to ROI manager and saved

run("Cell Counter");
waitForUser("Do this:  mark Keep original. Load caveolae markers. Run Measure. Close plugin");
selectWindow("Results");
saveAs("Results", dir+ title + "_" + cell+ "_PART_" + part+ "_POSICIONES_REANALYSIS" + ".txt");
run("Results... ", "open=["+ dir+ title+ "_" + cell+ "_PART_" + part+ "_POSICIONES_REANALYSIS.txt]");

selectWindow("Results");
numberOfPoints = getValue("results.count");
for (i = 0; i < numberOfPoints; i++) {
	XM = getResult("X", i);
	YM = getResult("Y", i);
	makePoint(XM, YM);
	roiManager("Add");
}


selectWindow("Results");
run("Close");
cav = roiManager("count"); 

//Manually open membranes
waitForUser("Open PM and LD membranes to ROI manager");
selectWindow("ROI Manager");
roiManager("Save", dir+  title+ "_" +cell+ "_PART_" + part+ "_CAVEOLAEDISTANCES_REANALYSIS.zip");

//Membranes are the last ROIs. We now obtain minimum distance from each caveolae to each membrane
//the sum of both distances will be the through caveolae PM-LD distance

run("Set Scale...", "distance=676.4536 known=1 pixel=1 unit=µm global");
getPixelSize(unit, pixelWidth, pixelHeight);
x1 = roiManager("Count");
//roiManager("Select", 0);
//current = Roi.getName();
line = x1;
s = 1; 
e = x1-2; //this is LD membrane
line2= x1-1; //this is PM


//calculate minimum distance caveolae-PM
roiManager("Select", line-1); //we select PM
run("Interpolate", "interval=1 smooth");
getSelectionCoordinates(linex, liney);
Array.show(linex, liney);
selectWindow("Arrays");
saveAs("Results", dir+title+ "_" + cell+ "_PART_"  +part+"_" + "PMCOORDINATES1_REANALYSIS"  + ".txt");
selectWindow(title+ "_" + cell+ "_PART_"  +part+"_" + "PMCOORDINATES1_REANALYSIS"  + ".txt");
run("Close");

DISTtoCM = newArray(e-(s-1));
ArraycoordminCPMX= newArray(e-(s-1));
ArraycoordminCPMY= newArray(e-(s-1));
ArraycoordminCLDX= newArray(e-(s-1));
ArraycoordminCLDY= newArray(e-(s-1));
CoordCX= newArray(e-(s-1));
CoordCY= newArray(e-(s-1));

for(i = s-1, t = 0; i < e; i++, t++){
roiManager("Select", i);
getSelectionCoordinates(xc, yc);
min = 100000;
for(j = 0; j < linex.length; j++){
for(k = 0; k < xc.length; k++){
dist = calculateDistance(linex[j], liney[j], xc[k], yc[k]);
if(dist < min){
min = dist;
coordminCPMX=linex[j];
coordminCPMY=liney[j];
}
								}
					}
makeLine(coordminCPMX,coordminCPMY, xc[0], yc[0]);
run("Add Selection...");
DISTtoCM[t] = min*pixelWidth;
ArraycoordminCPMX[t]= coordminCPMX;
ArraycoordminCPMY[t]= coordminCPMY;
CoordCX[t]= xc[0];
CoordCY[t]= yc[0];

}

//calculate minimum distance caveolae-LD
roiManager("Select", line2-1);
run("Interpolate", "interval=1 smooth");
getSelectionCoordinates(line2x, line2y);
Array.show(line2x, line2y);
selectWindow("Arrays");
saveAs("Results", dir+title+ "_" + cell+ "_PART_"  +part+"_" + "LDCOORDINATES1_REANALYSIS"  + ".txt");
selectWindow(title+ "_" + cell+ "_PART_"  +part+"_" + "LDCOORDINATES1_REANALYSIS"  + ".txt");
run("Close");

DISTtoLD = newArray(e-(s-1));
for(i = s-1, t = 0; i < e; i++, t++){
roiManager("Select", i);
getSelectionCoordinates(xc, yc);
min = 100000;
for(j = 0; j < line2x.length; j++){
for(k = 0; k < xc.length; k++){
dist = calculateDistance(line2x[j], line2y[j], xc[k], yc[k]);
if(dist < min){
min = dist;
coordminCLDX=line2x[j];
coordminCLDY=line2y[j];
}
								}
					}
makeLine(coordminCLDX,coordminCLDY, xc[0], yc[0]);
run("Add Selection...");
DISTtoLD[t] = min*pixelWidth;
ArraycoordminCLDX[t]= coordminCLDX;
ArraycoordminCLDY[t]= coordminCLDY;
}

roiManager("Deselect");
run("To ROI Manager");
roiManager("Save", dir + title+ "_" +cell+ "_PART_" + part+ "_ROISDISTANCESCAVEOLAE_REANALYSIS.zip");
selectWindow("ROI Manager");
run("Close");
roiManager("Open", dir + title+ "_" +cell+ "_PART_" + part+ "_CAVEOLAEDISTANCES_REANALYSIS.zip");

//We add caveolae coordinates to results

run("Set Measurements...", "  redirect=None decimal=9");
arrayRM = newArray("0");;
for (i=1;i<roiManager("count");i++){
arrayRM = Array.concat(arrayRM,i);
}
roiManager("select", arrayRM);
roiManager("Measure");
selectWindow("Results");
//print(lengthOf(DISTtoLD));

//We transform results columns to arrays
XCoord=newArray(x1-2);
for (i=0; i<lengthOf(DISTtoLD); i++){
	XCoord[i] = getResult("X",i);
	}

YCoord=newArray(x1-2);
for (i=0; i<lengthOf(DISTtoLD); i++){
	YCoord[i] = getResult("Y",i);
}
//Array.show(XCoord);
//Array.show(YCoord);


//distance results	
selection = Array.getSequence(DISTtoLD.length);
for(i = 0; i < selection.length; i++)
selection[i] = selection[i] + 1;

//all results together in table
Array.show(selection,DISTtoCM, DISTtoLD, XCoord, YCoord,ArraycoordminCPMX,ArraycoordminCPMY,ArraycoordminCLDX,ArraycoordminCLDY,CoordCX,CoordCY);


//save and close	
selectWindow("Arrays");
saveAs("Results", dir + title + "_" +cell+ "_PART_" + part  + "_CAVEOLAEDISTANCES_REANALYSIS" + ".txt");
roiManager("deselect");

selectWindow("Results");
run("Close");
selectWindow(title + "_" +cell+ "_PART_" + part+ "_CAVEOLAEDISTANCES_REANALYSIS" + ".txt");
run("Close");


//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////
///////PART1.2 REMOVE CAVEOLAE AND LEAVE ONLY MEMBRANES/////////////
////////////////////////////////////////////////////////////////////////////////////


rois = newArray(cav);
for (i=0; i<rois.length; i++)
      rois[i] = i;



roiManager("Select", rois);
roiManager("Delete");


//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////
///////PART 2 INTERPOLATE MEMBRANES INTERVAL = 8//////////////
////////////////////////////////////////////////////////

title = getTitle;
dir = getDirectory("image");
run("Set Scale...", "distance=676.4536 known=1 pixel=1 unit=µm");



getPixelSize(unit, pixelWidth, pixelHeight);
membrane = 1;
LD = 0;

roiManager("Select", membrane);
run("Interpolate", "interval=8");
roiManager("Add");
roiManager("Select", 2);
getSelectionCoordinates(linex, liney);
Array.show(linex, liney);
selectWindow("Arrays");
saveAs("Results", dir+title+ "_" + cell+ "_PART_"  +part+"_" + "PMCOORDINATES8_REANALYSIS"  + ".txt");
selectWindow(title+ "_" + cell+ "_PART_"  +part+"_" + "PMCOORDINATES8_REANALYSIS"  + ".txt");
run("Close");
	
roiManager("Select", LD);
run("Interpolate", "interval=8");
roiManager("Add");
roiManager("Select", 3);
getSelectionCoordinates(line2x, line2y);
Array.show(line2x, line2y);
selectWindow("Arrays");
saveAs("Results", dir+title+ "_" + cell+ "_PART_"  +part+"_" + "LDCOORDINATES8_REANALYSIS"  + ".txt");
selectWindow(title+ "_" + cell+ "_PART_"  +part+"_" + "LDCOORDINATES8_REANALYSIS"  + ".txt");
run("Close");


roiManager("Save", dir+title+  "_" + cell+ "_PART_" + part +"_MEMBRANEROIS8_REANALYSIS" + ".zip");

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////
///////PART 3 CALCULATE GLOBAL PM-LD DISTANCES//////////////
////////////////////////////////////////////////////////

//Get image name and directory
title = getTitle;
dir = getDirectory("image");
run("Set Scale...", "distance=676.4536 known=1 pixel=1 unit=µm");
getPixelSize(unit, pixelWidth, pixelHeight);

roiManager("Select", 2);
getSelectionCoordinates(linex, liney);
results = newArray(linex.length);
PMX= newArray(linex.length);
PMY= newArray(linex.length);
closerLDX= newArray(linex.length);
closerLDY= newArray(linex.length);

	
roiManager("Select", 3);
getSelectionCoordinates(xc, yc);


for(j = 0; j < linex.length; j++){
	min = 100000;
	for(k = 0; k < xc.length; k++){
		dist = calculateDistance(linex[j], liney[j], xc[k], yc[k]);
		if(dist < min){
			min = dist;
			coordmin1PMX=linex[j];
			coordmin1PMY=liney[j];
			coordmin1LDX=xc[k];
			coordmin1LDY=yc[k];
			}
	}
	makeLine(coordmin1PMX,coordmin1PMY, coordmin1LDX, coordmin1LDY);
	run("Add Selection...");
	results[j] = min;
	//results[j] = min*pixelWidth;
	PMX[j]=coordmin1PMX;
	PMY[j]=coordmin1PMY;
	closerLDX[j]=coordmin1LDX;
	closerLDY[j]=coordmin1LDY;
}
roiManager("Deselect");
run("To ROI Manager");
roiManager("Save", dir+ title+ "_" +cell+ "_PART_" + part+ "_ROISALLDISTANCESPMLD_REANALYSIS.zip");
selectWindow("ROI Manager");
run("Close");
roiManager("Open", dir+  title+ "_" +cell+ "_PART_" + part+ "_MEMBRANEROIS8_REANALYSIS.zip");
selectWindow(title);
Array.show(results,PMX,PMY, closerLDX,closerLDY);

//save and close
	selectWindow("Arrays");
	saveAs("Results", dir+title+ "_" + cell+ "_PART_" + part  +"_" + "ALLDISTANCESPMLD_REANALYSIS"  + ".txt");
	selectWindow(title+ "_" + cell+ "_PART_" + part  +"_" + "ALLDISTANCESPMLD_REANALYSIS"  + ".txt");
	run("Close");


results = newArray(xc.length);
closerPMX= newArray(xc.length);
closerPMY= newArray(xc.length);
LDX= newArray(xc.length);
LDY= newArray(xc.length);

for(j = 0; j < xc.length; j++){
	min = 100000;
	for(k = 0; k < linex.length; k++){
		dist = calculateDistance(xc[j], yc[j], linex[k], liney[k]);
		if(dist < min){
			min = dist;
			coordmin1PMX=linex[k];
			coordmin1PMY=liney[k];
			coordmin1LDX=xc[j];
			coordmin1LDY=yc[j];
			}
	}
	makeLine(coordmin1PMX,coordmin1PMY, coordmin1LDX, coordmin1LDY);
	run("Add Selection...");
	results[j] = min;
	//results[j] = min*pixelWidth;
	closerPMX[j]=coordmin1PMX;
	closerPMY[j]=coordmin1PMY;
	LDX[j]=coordmin1LDX;
	LDY[j]=coordmin1LDY;
}

roiManager("Deselect");
run("To ROI Manager");
roiManager("Save", dir+ title+ "_" +cell+ "_PART_" + part+ "_ROISALLDISTANCESLDPM_REANALYSIS.zip");
selectWindow("ROI Manager");
run("Close");
roiManager("Open", dir+  title+ "_" +cell+ "_PART_" + part+ "_MEMBRANEROIS8_REANALYSIS.zip");
selectWindow(title);
Array.show(results,LDX, LDY, closerPMX,closerPMY);

//save and close
	selectWindow("Arrays");
	saveAs("Results", dir+title+ "_" + cell+ "_PART_" + part  +"_" + "ALLDISTANCESLDPM_REANALYSIS"  + ".txt");
	selectWindow(title+ "_" + cell+ "_PART_" + part  +"_" + "ALLDISTANCESLDPM_REANALYSIS"  + ".txt");
	run("Close");


//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
///////PART 4: MEASURE PERIMETER//////////////
////////////////////////////////////////////////////////


run("Set Scale...", "distance=676.4536 known=1 pixel=1 unit=µm");
dir= getDirectory("image");
title=getTitle();
//Dialog.create("Imagen");
//Dialog.addString("Célula","");
//Dialog.show();
//cell=Dialog.getString;

roiManager("Select", 2);
roiManager("Measure");
selectWindow("Results");
saveAs("Results", dir+title+ "_" + cell +"_PART_" + part  + "_" +"PERIMETROFINAL_REANALYSIS.txt");
//waitForUser("Apunta");
selectWindow("Results");
run("Close");
//selectWindow("ROI Manager");
//run("Close");


//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////
//////PART 5: LINES TO POLYGON FOR CYTOPLASMIC COORDINATES////
////////////////////////////////////////////////////////
//Run these lines to start macro here omitting previous steps
//title = getTitle;
///dir = getDirectory("image");
//roiManager("Select", 0);
//current = Roi.getName();
//Dialog.create("Imagen");
//Dialog.addString("Célula","");
//Dialog.show();
//cell=Dialog.getString;

//Dialog.create("Imagen");
//Dialog.addString("Parte","");
//Dialog.show();
//part=Dialog.getString;


//Open PM-LD membranes before running
//ROI 2 es PM. ROI 3 es LD

//roiManager("Select", 0);
//run("Interpolate", "interval=1 smooth");
//roiManager("Add");
//roiManager("Select", 1);
//run("Interpolate", "interval=1 smooth");
//roiManager("Add");


///CHECK FIRST IF THE POLYGON IS OK
//delete interpolated membranes
ultimaroi = roiManager("Count");
roiManager("Select", ultimaroi-1);
roiManager("Delete");
ultimaroi = roiManager("Count");
roiManager("Select", ultimaroi-1);
roiManager("Delete");

roiManager("Select", 1);
run("Interpolate", "interval=1 smooth");
roiManager("Add");
roiManager("Select", 0);
run("Interpolate", "interval=1 smooth");

roiManager("Add");
roiManager("Select", 3);
getSelectionCoordinates(ldx, ldy);
roiManager("Select", 2);
//Array.show(ldx);
//Array.show(ldy);
getSelectionCoordinates(mpx, mpy);
//coordx=Array.concat(mpx,Array.reverse(ldx)) ;
//coordy=Array.concat(mpy,Array.reverse(ldy)) ;

coordx=Array.concat(mpx,ldx) ;
coordy=Array.concat(mpy,ldy) ;

//Array.show(coordx, coordy);
makeSelection( "polygon", coordx, coordy );

//check that the polygon is ok, if not, it will do alternative
label="Polygon OK?";
Dialog.create("Polygon");
Dialog.addCheckbox(label, true);
Dialog.show();
answer=Dialog.getCheckbox;


if (answer==1) {

roiManager("Add");
roiManager("Save", dir + title+ "_" +cell+ "_PART_" + part+ "_POLYGON_REANALYSIS.zip");


} else {

roiManager("Select", 3);
getSelectionCoordinates(ldx, ldy);
roiManager("Select", 2);
//Array.show(ldx);
//Array.show(ldy);
getSelectionCoordinates(mpx, mpy);
coordx=Array.concat(mpx,Array.reverse(ldx)) ;
coordy=Array.concat(mpy,Array.reverse(ldy)) ;

//coordx=Array.concat(mpx,ldx) ;
//coordy=Array.concat(mpy,ldy) ;

//Array.show(coordx, coordy);
makeSelection( "polygon", coordx, coordy );
roiManager("Add");
roiManager("Save", dir + title+ "_" +cell+ "_PART_" + part+ "_POLYGON_REANALYSIS.zip");
	
}




////////////////////////////////////////////////////////
//////PART 6: All pixel coordinates////
////////////////////////////////////////////////////////

title = getTitle;
dir = getDirectory("image");

getDimensions(width, height, channels, slices, frames);
newImage(title+ "_" +"Mask" + "_" +cell + "_"+part, "8-bit white", width, height, slices);
title2 = getTitle;

poligono = roiManager("Count")-1;
roiManager("Select", poligono);
run("Fill", "slice");
run("Save XY Coordinates...", "background=255 invert save=" + dir + title + "_" + cell+ "_PART_" + part+  "_PIXELCOORDINATES_REANALYSIS.txt");
run("Results... ", "open="+ dir + title + "_" + cell+ "_PART_" + part+ "_PIXELCOORDINATES_REANALYSIS.txt");
selectWindow("Results");
saveAs("Results", dir+title+ "_" + cell+ "_PART_"  +part + "_PIXELCOORDINATES_REANALYSIS"  + ".txt");
selectWindow("Results");
run("Close");
selectWindow(title2);
run("Close");
selectWindow("Log");
run("Close");
selectWindow("ROI Manager");
run("Close");
//run("Close All");

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//FUNCTIONS

function calculateDistance(x1, y1, x2, y2){
	temp = sqrt(pow((x2-x1), 2) + pow((y2-y1),2));
	return temp;
}

