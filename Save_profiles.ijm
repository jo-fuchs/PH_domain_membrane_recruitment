///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Membrane recruitment analysis part 1
//
// Run for every image separately after:
//	 - opening Image-Stack with 4 channels
//	 - Lines stored in ROI-manager
//
// Goals:
//   - check &/or create LineProfiles-folder
//   - Create line profile of a stored line in the ROI manager in all 4 channels
//   - Scale Y-axis to µm instead of pixeles
//	 - Save final table for each ROI with the name: "ImageName_RoiName"
//	 - Save ROIs
//
//   Version 0.1 (09.04.2020)
// 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// select working directory & read image name & ID
path = getDirectory("Choose a Directory");
imgID = getImageID;
imgName = getTitle();;

// create Results directory
subdir=path+"LineProfiles/";
File.makeDirectory(subdir); 


// loop through Roi-Manager for each cell
nR = roiManager("Count"); 

for (k=0; k<nR; k++) { 
	run("Clear Results");
	roiManager("Select",k);
	rName = Roi.getName;
	Stack.setChannel(2);
	profile_Mem = getProfile();		// get profile for membrane marker
	Stack.setChannel(1);
	profile_PH = getProfile();		// get profile for Domain of interest
	
	Stack.setChannel(3);
	profile_Dapi = getProfile();
	Stack.setChannel(4);
	profile_pAKT = getProfile();

	// scale X-axis from pixel to scale
	x = Array.getSequence(profile_Mem.length);
	scaleArray(imgID, x );

	//create results table with  X and Y of Profiles
	//// if you want to give different names make sure to also modify them accordingly in the R-code (part 2) of the analysis!)
	for (i=0; i<profile_Mem.length; i++) {
		setResult("X", i, x[i]);
		setResult("Actin", i, profile_Mem[i]) ;
		setResult("PH", i, profile_PH[i]) ;
		setResult("Dapi", i, profile_Dapi[i]) ;
		setResult("pAKT", i, profile_pAKT[i]) ;
	}
	// save results table as "ImageName_RoiName.txt"
	updateResults(); 
	saveAs("Measurements", subdir+imgName+"_"+rName+".txt"); 
	}
	
	// Save ROIs as "ImageName.zip"
     roiManager("deselect");
     roiManager("save", path+imgName+".zip") 



//// scale Array function from user anon96376101 at https://forum.image.sc/t/getprofile-function-not-following-the-scale/9260
function scaleArray( image, array ) {
	front = getImageID;
	selectImage( image );
	for ( i= 0; i < array.length; i++ ) toScaled( array[i] );
	selectImage( front );
}

