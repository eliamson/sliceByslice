macro "sliceBySlice" {

//V0.9.4 rename
//Eli Amson eli.amson1988@gmail.com

//This macro measures slice by slice area and compactness of object of interest. It should be applied to a stack (grey-values or already binarized) of which the unit is either mm or inch.

run("Set Measurements...", "area area_fraction redirect=None decimal=3");

//Uses "Optimise Threshold" of BoneJ

	//If stack not already in mm, converts inch to mm (error if neither)
	getPixelSize(unit, pw, ph, pd);
	if (unit=="inches"){
	setVoxelSize(pw*25.4, ph*25.4, pd*25.4,"mm");
	getPixelSize(unit, pw, ph, pd);
	}else{
		if(unit!="mm"){
			exit("Not in inches or mm");
		}
	}

	// Binarize stack if not already done
	if (!is("binary")){
	run("Optimise Threshold", "threshold apply");
	}else{
		if (getPixel(0, 0)==255){//In case inverted binarization
			run("Invert", "stack");
			waitForUser("Top left pixel was not background, so stack was inverted");
			}
	}
	// Close previous Results window if any
	listWin = getList("window.titles");
	  if (listWin.length!=0){
		for (i=0;i<listWin.length;i++){
			if (listWin[i]=="Results"){
			selectWindow("Results");
			run("Close");
				} 
			}
	  }
	// If a minimum  particle size is to be set (here chosen as 10% of width of pic). 
	getPixelSize(unit, pw, ph, pd);
	minPart=pow(getWidth()*pw*0.1,2);

	// Area and compactness of largest particle (i.e., studied bone, larger than min size, see above) of each slice 
	ResArea=newArray(nSlices);
	ResC=newArray(nSlices);

	for (i=1;i<=nSlices;i++){
		setSlice(i);
		run("Analyze Particles...", "size=0-Infinity show=Nothing clear  include slice");
		maxPartN=0;
		maxPartArea=0;
			if (nResults!=0){	
				for (j=0;j<nResults;j++){
					if (getResult("Area",j)>maxPartArea){
					maxPartN=j;
					maxPartArea=getResult("Area",j);
					}
				}
			ResArea[i-1]=getResult("Area",maxPartN);
			ResC[i-1]=getResult("%Area",maxPartN);
			}else{
			ResArea[i-1]=NaN;
			ResC[i-1]=NaN;
			}
			run("Clear Results");
	}
	
	// Make final Result table
	for (i=0; i<nSlices; i++) {
		setResult("ResArea ("+unit+")", i, ResArea[i]); // ResArea = whole sectional area of the studied bone (on each slice)
		setResult("ResC", i, ResC[i]) ; // ResC = global compactness of the studied bone's cross-section (on each slice)
	}
	// Selection of first (Z1) and last (Z2) slices and slices to remove (semi-automatically)
	setSlice(1);
	setTool("wand");
	run("Wand Tool...", "tolerance=0 mode=8-connected");
	waitForUser("Move to first slice with complete bone (Z1)");
	Z1=getSliceNumber();
	ToDel=true;
	ToRemCount=0;
		while (ToDel){
		waitForUser("Check if slices to remove and then click 'OK'");
		Dialog.create("Slices to exclude?");
		Dialog.addCheckbox("Are there (more) slices to exclude?", false)
		Dialog.show()
		ToDel=Dialog.getCheckbox();
			if (ToDel) {
			waitForUser("Move to first slice to remove");
			F=getSliceNumber();
			waitForUser("Move to last slice to remove");
			L=getSliceNumber();
			print("Last slice excluded: "+L);
			if (F==L) {
				setResult("ToRem", L-1, 1) ;
				ToRemCount=ToRemCount+1;
				}else{
					for (j=F-1; j<L; j++) {
						ToRemCount=ToRemCount+1;
						setResult("ToRem", j, 1) ;
						}
					}
			}	
		}
	waitForUser("Move to last slice with complete bone (Z2)");
	Z2=getSliceNumber();
	setResult("ToRem", Z1-1,"Z1") ;
	setResult("ToRem", Z2-1,"Z2") ;

	//Print overall stats
	ResCmean=0;
	ResAreaMean=0;
		for (i=Z1-1; i<Z2; i++){	//Take range within Z1 and Z2 and exclude ToRems
			if (getResult("ToRem",i)!=1) {
				ResCmean=ResCmean+ResC[i];
				ResAreaMean=ResAreaMean+ResArea[i];
				}
			}
		Nincluded=Z2-Z1+1-ToRemCount;
		ResCmean=ResCmean/Nincluded;
		ResAreaMean=ResAreaMean/Nincluded;
		print("Mean total cross-sectional area (Tt.Ar): "+ResAreaMean+" mm2");
		print("Mean global compactness(Cg): "+ResCmean+"%");
}

//ResAreaCl=newArray(Z2-Z1+1-ToRemCount);
//ResCCl=newArray(Z2-Z1+1-ToRemCount);
