close("*");

//Ask user to choose the input and output directories
directory = getDirectory("Choose input directory");
fileList = getFileList(directory);
outputDirectory = getDirectory("Choose output directory");
run("Bio-Formats Macro Extensions");

manualSettings = getBoolean("Do you want to manually input the set points?");

for (i=0; i<fileList.length; i++) {
	file = directory + fileList[i];
	Ext.setId(file);

	//Measure number of series in file
	Ext.getSeriesCount(nStacks);
	
	//Counts number of channels in series
	Ext.getSizeC(sizeC);

	//Create a new array for saving the sample names, GFP bundle count and GFP bundle density of normalized stacks
	sampleNameArray = newArray(nStacks);
	bundleRatioArrray = newArray(nStacks);
	bundleCountArray = newArray(nStacks);

	//Initialize variables that need to be checked on the first iteration
	if(manualSettings || i == 0){
		lowestChannel = -1;
		lowerThreshold = -1;
	}
	
	//Open all stacks from set of lif files, one stack at a time
	for(j=0; j<nStacks; j++) {	

		//Once the first run is complete turn on batch mode
		if(lowestChannel == -1 && manualSettings){
			setBatchMode(false);
		}
		else{
			setBatchMode(true);
		}

		//Show/update progress to user in a bar 
		progress = (j*(i+1)-1)/(fileList.length*nStacks);
		showProgress(progress);

		//Only open the image if it has one / (i.e. part of collection, but not subcollection), and ends with a series number (not a, b, merged, etc.)
		Ext.setSeries(j);
		Ext.getSeriesName(seriesTitle)
		print(seriesTitle);
		if(!matches(seriesTitle, ".*/.*/.*$") && matches(seriesTitle, ".*Series[0-9][0-9][0-9]$")){
			run("Bio-Formats Importer", "open=file color_mode=Default split_channels view=[Standard ImageJ] stack_order=Default series_"+d2s(j+1,0)); 
	
			imageTitle = getTitle();
			print(imageTitle);
	
			//if this is the first stack in a lif file, find dimmest channel and ask user if this is the correct 
			if(lowestChannel == -1){
				for(a=0; a<sizeC; a++){
					currentTitle = replace(imageTitle, "C=.", "C=" + a);
					selectWindow(currentTitle);
					Stack.getStatistics(dummy, mean);
					if(a == 0){
						lowestMean = mean;
						lowestChannel = a;
					}
					if(mean < lowestMean){
						lowestMean = mean;
						lowestChannel = a;
					}
				}
				currentTitle = replace(imageTitle, "C=.", "C=" + lowestChannel);
				selectWindow(currentTitle);
				setSlice(round(nSlices/2));
				if(manualSettings){
					run("Enhance Contrast", "saturated=0.35");
					GFPchannel = getNumber("Which channel contains the volume objects you want to measure?", lowestChannel);
				}
			}
	
			//Close all channels that are not the desired channel
			for(a=0; a<sizeC; a++){
				currentTitle = replace(imageTitle, "C=.", "C=" + a);
				selectWindow(currentTitle);
				if(a != lowestChannel){
					close();
				}
			}
	
			GFPstack = getTitle();
			//Set 3D object counter measurement parameters
			run("3D OC Options", "volume surface nb_of_obj._voxels nb_of_surf._voxels integrated_density mean_gray_value std_dev_gray_value median_gray_value minimum_gray_value maximum_gray_value centroid mean_distance_to_surface std_dev_distance_to_surface median_distance_to_surface centre_of_mass bounding_box dots_size=5 font_size=10 redirect_to=none");
			
			//Clean-up image and find optimal image threshold, then ask user to verify value
			run("Median 3D...", "x=3 y=3 z=3");
	
			//For the first iteration, find the optimal minimum cutoff to hold fixed for all subsequent analyses
			if(lowerThreshold == -1){
				setAutoThreshold("Triangle dark stack");
				getThreshold(lowerThreshold, dummy);
				if(manualSettings){
					lowerThreshold = getNumber("What is the minimum intensity cutoff you would like to use?", lowerThreshold);
				}
			}
			
			//Measure the GFP objects in the image
			run("3D Objects Counter", "threshold=" + lowerThreshold + " slice=1 min.=4 max.=93126656 exclude_objects_on_edges objects statistics summary");
			
			//Save the objects map
			selectWindow("Objects map of " + GFPstack);
			objectMapName = "Objects map of " + replace(GFPstack, fileList[i] + " - ", "");
			objectMapName = replace(objectMapName, "/Series... - C=.", " - " + j + ".tif");
			saveAs("tiff", outputDirectory + objectMapName);
			close();
			
			//Since the image stack started above the GFP and finished below the GFP, and the images are not level (one side higher than the other),
			//to be able to estimate the GFP density, the stack needs to be cropped to a normalized volume containing GFP
			//The definition used below is that slice on either end of the stack will be removed until a slice is reached that is >= the mean of the slices.
			
			//Measure the mean intensity of each slice within the stack, and save each measurement in a corresponding array
			selectWindow(GFPstack);
			sliceMeanArray = newArray(nSlices);
			for(a=1; a<=nSlices; a++){
				setSlice(a);
				getStatistics(dummy, sliceMean);
				sliceMeanArray[a-1] = sliceMean;
			}
			
			//Find the mean of the array to use as a reference when cropping the stack
			Array.getStatistics(sliceMeanArray, dummy, dummy, sliceArrayMean, dummy);
			
			//Look for the slice closest to the top and the slice closest to the bottom that have and intensity >= the mean
			for(a=0; a<sliceMeanArray.length; a++){
				if(sliceMeanArray[a] >= sliceArrayMean){
					meanTopSlice = a+1;
					a = sliceMeanArray.length;
				}
			}
			
			for(a=sliceMeanArray.length-1; a>=0; a--){
				if(sliceMeanArray[a] >= sliceArrayMean){
					meanBottomSlice = a+1;
					a = -1;
				}
			}
			
			//Rescale the upper and lower normalized stack limits to microns, and then search for all objects whose centroids are within these limits
			getVoxelSize(voxelWidth, voxelHeight, voxelDepth, unit);
			getDimensions(stackWidth, stackHeight, dummy, dummy, dummy);
			normalizeStackTop = meanTopSlice * voxelDepth;
			normalizeStackBottom = meanBottomSlice * voxelDepth;
			
			selectWindow("Results");
			normalizedTotalCount = 0;
			normalizedTotalVolume = 0;
	
			for(a=0; a<nResults; a++){
				gfpZposition = getResult("ZM", a);
				gfpVolume = getResult("Volume (micron^3)", a);
				if(gfpZposition >= normalizeStackTop && gfpZposition <= normalizeStackBottom){
					normalizedTotalCount++;
					normalizedTotalVolume = normalizedTotalVolume + gfpVolume;
				}
			}
	
			//Calculate the density and count of the samples 
			normalizedStackVolume = (stackWidth * voxelWidth * stackHeight * voxelHeight * (normalizeStackBottom - normalizeStackTop));
			normalizedPercentVolume = normalizedTotalVolume / normalizedStackVolume;
			countPerMicroLiter = (normalizedTotalCount * 1000000000 / normalizedStackVolume);
	
			//Store the results in the corresponding array
			bundleRatioArrray[j] = normalizedPercentVolume;
			bundleCountArray[j] = countPerMicroLiter;
	
			//Store the sample name in the sample name array
			sampleName = replace(GFPstack, fileList[i] + " - ", "");
			sampleNameArray[j] = replace(sampleName, "/Series... - C=.", "");
			
			//Save the results table
			selectWindow("Results");
			statisticsName = "Statistics for " + replace(GFPstack, fileList[i] + " - ", "");
			statisticsName = replace(statisticsName, "/Series... - C=.", " - " + j + ".xls");
			saveAs("Results", outputDirectory + statisticsName);
			run("Close");
			close("*");
		}

	}
	
	//Print the density and count statistics into a speadsheet
	for(a=0; a<nStacks; a++){
		setResult("sampleID", a, sampleNameArray[a]);
		setResult("GFP Volume Ratio", a, bundleRatioArrray[a]);
		setResult("GFP Count per uL", a, bundleCountArray[a]);
	}

	//Save the resulting table
	tableName = "Normalized density statistics for " + replace(fileList[i], ".lif", ".xls");
	saveAs("Results", outputDirectory + tableName);
	run("Close");
}

setBatchMode(false);
