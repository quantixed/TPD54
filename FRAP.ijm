/*
 * Do simple FRAP analysis
 * Ultraview hyperstack with data in channel 1 and single FRAP ROI in channel 2
 * For 0.069 um/pixel
 * The code makes an ROI for the FRAP region (specified in the file as half-sze)
 * It asks the user to make a bg selection and then specify the cell.
 * The csv output is saved to the directory the user choses, using the name of the file.
 * The ROIs are saved as a zip for reproducibility
 * Written by Stephen Royle 2018-11-17
 */
function lowestVal(colStr)	{
	theLowest = 10000	// set to something big
	for (row = 0; row < nResults; row++) { 
        theLowest = minOf(theLowest,getResult(colStr,row));
	}
	return theLowest
}
pxSize = 0.068964788;
// get the name of the original window
win = getTitle();
run("Set Measurements...", "area mean min bounding integrated display redirect=None decimal=3");
// make window for FRAP ROI
run("Duplicate...", "title=tempImg duplicate channels=2 frames=1");
// make data stack (remove FRAP channel)
selectWindow(win);
run("Duplicate...", "title=dataImg duplicate channels=1");
selectWindow(win);
close();
// find the FRAP ROI
selectWindow("tempImg");
setThreshold(0, 1);
setOption("BlackBackground", false);
run("Convert to Mask");
run("Analyze Particles...", "display clear summarize add");
// because the 0th row does not necessarily contain the top left ROI
bx=2 * round(lowestVal("BX") / pxSize);
by=2 * round(lowestVal("BY") / pxSize);
ww=2 * round(getResult("Width",0) / pxSize);
hh=2 * round(getResult("Height",0) / pxSize);
//print(bx,by,ww,hh);
roiManager("reset");
close("tempImg");
// now select the data window and make all the selections
selectWindow("dataImg");
run("Specify...", "width="+ww+" height="+hh+" x="+bx+" y="+by+" slice=1");
roiManager("Add");
// Now user defines the background
	setTool(0);
	waitForUser("Define background", "Make a box away from the cell");
roiManager("Add");
// Now user defines the background
	setTool(3);
	waitForUser("Define cell", "Draw around the cell");
roiManager("Add");
run("Clear Results");
roiManager("Multi Measure");
dir=getDirectory("Select a directory");
saveAs("Results", dir+win+".csv");
roiManager("Save", dir+win+".zip");