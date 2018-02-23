#pragma TextEncoding = "MacRoman"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

// Menu item for easy execution
Menu "Macros"
	"IMOD Models...",  IMODModelAnalysis()
	"Start Over", CleanSlate()
End

Function IMODModelAnalysis()
	LoadIMODModels()
	SetUpWindows()
	ProcessAllModels()
	CollectAllMeasurements()
//	MakeSummaryLayout()
//	RotatePits()
//	OverlayAllPits()
	FormatPlots()
	MakeTheLayouts()
End

Function LoadIMODModels()
	// Check we have FileName wave and PixelSize
	Wave/T/Z FileName = root:FileName
	Wave/Z PixelSize = root:PixelSize
	if (!waveexists(FileName))
		Abort "Missing FileName textwave"
	endif
	if(!WaveExists(PixelSize))
		Abort "Missing PixelWave numeric wave"
	endif
	
	NewDataFolder/O/S root:data
	
	String expDiskFolderName, expDataFolderName
	String FileList, ThisFile, pdfName
	Variable FileLoop, nWaves, i
	
	NewPath/O/Q/M="Please find disk folder" ExpDiskFolder
	if (V_flag!=0)
		DoAlert 0, "Disk folder error"
		Return -1
	endif
	PathInfo /S ExpDiskFolder
	ExpDiskFolderName=S_path
	FileList=IndexedFile(expDiskFolder,-1,".txt")
	Variable nFiles=ItemsInList(FileList)
	
	for (FileLoop = 0; FileLoop < nFiles; FileLoop += 1)
		ThisFile = StringFromList(FileLoop, FileList)
		expDataFolderName = ReplaceString(".txt",ThisFile,"")
		NewDataFolder/O/S $expDataFolderName
		LoadWave/A/J/D/O/K=1/V={" "," $",0,0}/L={0,0,0,1,0}/P=expDiskFolder ThisFile
		MakeObjectContourWaves()
		SetDataFolder root:data:
	endfor
End

Function MakeObjectContourWaves()
	Concatenate/O/KILL wavelist("wave*",";",""), matA
	WaveStats/Q/RMD=[][0] matA
	// Scale the coordinates to real values
	ScaleCoords(matA)
	Variable nObjects = V_max + 1
	Variable nContours, contourVar
	String wName
	
	Variable i,j
	
	for (i = 0; i < nObjects; i += 1)
		MatrixOP/O filtObj = col(matA,0)
		filtObj[] = (filtObj[p] == i) ? matA[p][1] : NaN
		WaveTransform zapnans filtObj
		FindDuplicates/RN=uniqueContours filtObj
		nContours = Numpnts(uniqueContours)
		// zero-indexed list of contours in this object
		for (j = 0; j < nContours; j += 1)
			contourVar = uniqueContours[j]
			// find the rows that correspond to each contour
			Duplicate/O/FREE matA,matB
			matB[][] = (matB[p][0] == i && matB[p][1] == contourVar) ? matB[p][q] : NaN
			MatrixOp/O/FREE xW = col(matB,2)
			MatrixOp/O/FREE yW = col(matB,3)
			// no need to take the z column here
			WaveTransform zapnans xW
			WaveTransform zapnans yW
			// Now make ObjectContour waves
			// Object 0 is Mitochondria, Object 1 is vesicles
			if (i == 0)
				wName = "Mt" // Mitochondria
			elseif (i == 1)
				wName = "Vs" // Vesicle
			else
				wName = "Uk" + num2str(i) // unknown
			endif
			wName += "_" + num2str(contourVar)
			Concatenate/O/NP=1 {xW,yW}, $wName
		endfor
	endfor
	KillWaves/Z matA,filtObj,UniqueContours
End

///	@param	matA	wave reference to matrix
Function ScaleCoords(matA)
	Wave matA
	
	String txtName = ReplaceString("'",GetDataFolder(0),"")
	Wave/T/Z FileName = root:FileName
	Wave/Z PixelSize = root:PixelSize
	Wave/Z matA
	Variable pxSize
	
	if (!WaveExists(FileName) || !WaveExists(PixelSize))
		DoAlert 0, "Cannot scale"
		Return -1
	endif
	FindValue/TEXT=txtName FileName
	if (V_Value == -1)
		Print txtName, "didn't scale"
	endif
	
	pxSize = PixelSize[V_Value]
	matA[][2,4] *= pxSize
	// Print txtName, pxSize
End

Function SetUpWindows()
	String windowList = "allVsPlot;allRotVsPlot;"
	Variable nWindows = ItemsInList(windowList)
	String PlotName
	
	Variable i
	
	for(i = 0; i < nWindows; i += 1)
		plotName = StringFromList(i,windowList)
		KillWindow/Z $plotName
		Display/N=$plotName
	endfor
End

Function ProcessAllModels()
	SetDataFolder root:data:	// relies on earlier load
	DFREF dfr = GetDataFolderDFR()
	String folderName
	Variable numDataFolders = CountObjectsDFR(dfr, 4)
	
	Variable i
		
	for(i = 0; i < numDataFolders; i += 1)
		folderName = GetIndexedObjNameDFR(dfr, 4, i)
		SetDataFolder ":'" + folderName + "':"
		// run functions from here for everything we want to do
		FindVesicleCentresAndPlotOut()
		MitoPerimAndArea()
		LookAtModels()
//		Distance2PM()
		SetDataFolder root:data:
	endfor
	SetDataFolder root:
End

Function FindVesicleCentresAndPlotOut()
	// finding the centre of each vesicle and placing coords into a wave called VsWave
	String wList = WaveList("Vs_*",";","")
	Variable nWaves = ItemsInList(wList)
	Variable nCol = 3 // x y and z
	Make/O/N=(nWaves,nCol) Img_VsCentre
	Make/O/N=(nWaves) Img_VsMinAxis,Img_VsMajAxis,Img_VsPerimeter,Img_VsArea
	String currentDF = GetDataFolder(0)
	String wName,tName
	
	Variable i,j
	
	for (i = 0; i < nWaves; i += 1)
		wName = StringFromList(i, wList)
		Wave w0 = $wName
		nCol = dimsize(w0,1)
		for (j = 0; j < nCol; j += 1)
			WaveStats/Q/M=1/RMD=[][j] w0
			Img_VsCentre[i][j] = V_avg
		endfor
		// make a nice tracename
		TName = CleanupName(currentDF,0) + "_" + wName
		AppendToGraph/W=allVsPlot w0[][1]/TN=$tName vs w0[][0]
		// offset to origin
		ModifyGraph/W=allVsPlot offset($tName)={-Img_VsCentre[i][0],-Img_VsCentre[i][1]}
		// find eigenvectors and rotate vesicle coords (also offset them)
		Wave w1 = FindEV(w0)
		TName = CleanupName(currentDF,0) + "_" + NameOfWave(w1)
		AppendToGraph/W=allRotVsPlot w1[][1]/TN=$tName vs w1[][0]
		Img_VsMinAxis[i] = VesicleAxisLength(w1,1)
		Img_VsMajAxis[i] = VesicleAxisLength(w1,0)
		Img_VsPerimeter[i] = FindLengthOfXYCoords(w1)
		MatrixOp/O/FREE w1c0 = col(w1,0)
		MatrixOp/O/FREE w1c1 = col(w1,1)
		Img_VsArea[i] = PolygonArea(w1c0,w1c1)
	endfor
	if(numpnts(Img_VsArea) > 1)
		MatrixOp/O Img_VsAspectRatio = Img_VsMinAxis / Img_VsMajAxis
		MatrixOp/O Img_VsCircularity = (4 * pi * Img_VsArea) / (Img_VsPerimeter * Img_VsPerimeter)
	elseif(numpnts(Img_VsArea) == 1)
		MatrixOp/O Img_VsAspectRatio = Img_VsMinAxis / Img_VsMajAxis
		MatrixOp/O/FREE tempMat = Img_VsPerimeter * Img_VsPerimeter
		MatrixOp/O Img_VsCircularity = 4 * pi * tempMat
	else
		Print "No vesicles in", currentDF
	endif
End

///	@param	m1	2D wave of xy coords
///	@param	colNo	column number to use for search
Function VesicleAxisLength(m1,colNo)
	Wave m1
	Variable colNo
	MatrixOp/O/FREE m1c0 = col(m1,0)
	MatrixOp/O/FREE m1c1 = col(m1,1)
	Variable V_Value,len
	if(colNo == 1)
		FindLevel/Q/EDGE=1/P m1c0, 0
		len = abs(m1c1(V_Value))
		FindLevel/Q/EDGE=2/P m1c0, 0
		len += abs(m1c1(V_Value))
	else
		FindLevel/Q/EDGE=1/P m1c1, 0
		len = abs(m1c0(V_Value))
		FindLevel/Q/EDGE=2/P m1c1, 0
		len += abs(m1c0(V_Value))
	endif
	return len
End

///	@param	m1	2D wave of xy coords
Function FindLengthOfXYCoords(m1)
	Wave m1
	// make new 2D wave of xy coords
	Duplicate/O/FREE m1,tempDist
	// offset to zero
	tempDist[][0] -= m1[0][0]
	tempDist[][1] -= m1[0][1]
	// Differentiate, backward difference
	Differentiate/METH=2 tempDist
	// find norm, cumlative distance
	MatrixOp/O/FREE tempNorm = sqrt(sumRows(tempDist * tempDist))
	tempNorm[0] = 0 // first point is garbage
	// return the sum of distances
	return sum(tempNorm)
End

///	@param	m1	2D wave of xy coords
Function/WAVE FindEV(m1)
	Wave m1
	MatrixOp/O xCoord = col(m1,0)
	MatrixOp/O yCoord = col(m1,1)
	
	// translate to origin
	Variable offX = mean(xCoord)
	Variable offY = mean(yCoord)
	xCoord[] -= offX
	yCoord[] -= offY
	// do PCA. Rotated points are in M_R
	PCA/ALL/SEVC/SRMT/SCMT xCoord,yCoord
	WAVE M_R
	String mName = NameOfWave(m1) + "_r"
	Duplicate/O M_R, $mName
	Wave m2 = $mName
	// now thread it so the segment is contiguous
	InsertPoints DimSize(m2,0), 1, m2
	m2[DimSize(m2,0)-1][] = m2[0][q]
	Return m2
End

Function MitoPerimAndArea()
	// finding the centre of each vesicle and placing coords into a wave called VsWave
	String wList = WaveList("Mt_*",";","")
	Variable nWaves = ItemsInList(wList)
	Make/O/N=(nWaves) Img_MtPerimeter,Img_MtArea
	String currentDF = GetDataFolder(0)
	String wName,tName
	
	Variable i
	
	for (i = 0; i < nWaves; i += 1)
		wName = StringFromList(i, wList)
		Wave w0 = $wName
		Duplicate/O/FREE w0, tempMat
		InsertPoints DimSize(tempMat,0), 1, tempMat
		tempMat[DimSize(tempMat,0)-1][] = tempMat[0][q]
		Img_MtPerimeter[i] = FindLengthOfXYCoords(tempMat)
		MatrixOp/O/FREE w1c0 = col(tempMat,0)
		MatrixOp/O/FREE w1c1 = col(tempMat,1)
		Img_MtArea[i] = PolygonArea(w1c0,w1c1)
	endfor
End

Function FormatPlots()
	String windowList = "allVsPlot;allRotVsPlot;"
	Variable nWindows = ItemsInList(windowList)
	String PlotName
	
	Variable i
	
	for(i = 0; i < nWindows; i += 1)
		plotName = StringFromList(i,windowList)
		ModifyGraph/W=$plotName rgb=(29524,1,58982,6554)
		SetAxis/W=$plotName left -50,50
		SetAxis/W=$plotName bottom -50,50
		ModifyGraph/W=$plotName width={Aspect,1}
		ModifyGraph/W=$plotName mirror=1
	endfor
End

Function LookAtModels()
	// Plotting function to see what we have
	String wList = WaveList("Vs_*",";","")
	Variable nWaves = ItemsInList(wList)
	String currentDF = GetDataFolder(0)
	String plotName = "mod_" + CleanupName(currentDF,0)
	KillWindow/Z $plotName
	Display/N=$plotName/HIDE=1
	String wName,tName
	
	Variable i
	
	for (i = 0; i < nWaves; i += 1)
		wName = StringFromList(i, wList)
		if(Stringmatch(wName,"*_r") == 0) // filter out rotated objects
			Wave w0 = $wName
			// make a nice tracename - not really necessary
			tName = CleanupName(currentDF,0) + "_" + wName
			AppendToGraph/W=$plotName w0[][1]/TN=$tName vs w0[][0]
			ModifyGraph/W=$plotName rgb($tName)=(65535,0,65535) // vesicles are purple
		endif
	endfor
	wList = WaveList("Mt_*",";","")
	nWaves = ItemsInList(wList)
	for (i = 0; i < nWaves; i += 1)
		wName = StringFromList(i, wList)
		Wave w0 = $wName
		// make a nice tracename
		tName = CleanupName(currentDF,0) + "_" + wName
		AppendToGraph/W=$plotName w0[][1]/TN=$tName vs w0[][0]
		ModifyGraph/W=$plotName rgb($tName)=(32767,32767,32767) // mitos are grey
	endfor
	// Images are 1376 x 1032
	SetAxis/W=$plotName bottom 0,2000
	SetAxis/W=$plotName left 2000,0
	// scaling will be variable
	ModifyGraph/W=$plotName width={Plan,1,bottom,left}
	ModifyGraph/W=$plotName tick=3,mirror=1,noLabel=2,standoff=0
	ModifyGraph/W=$plotName margin=6
End

Function Distance2PM()
	Wave/Z VsWave
	String wList = WaveList("Mt_*",";","")
	String wName = StringFromList(0,wList) // it is probably PM_0
	Wave PMw = $wName
	
	Variable nFT = DimSize(VsWave,0)
	Make/O/N=(nFT) distWave,rowPM
	Variable i
	
	for (i = 0; i < nFT; i += 1)
		Duplicate/O PMw, m0
		m0[][] -= VsWave[i][q]
		MatrixOP/O result = m0 * m0
		MatrixOP/O result2 = sumrows(result)
		MatrixOP/O result3 = sqrt(result2)
		WaveStats/Q result3
		distWave[i] = V_min
		rowPM[i] = V_minLoc
	endfor
	KillWaves m0,result,result2,result3
End

Function CollectAllMeasurements()
	SetDataFolder root:data:	// relies on earlier load
	DFREF dfr = GetDataFolderDFR()
	String folderName
	Variable numDataFolders = CountObjectsDFR(dfr, 4)
	String wList = ""
	
	Variable i,j
	// assemble a string of semi-colon separated targets in the data folder
	for(i = 0; i < numDataFolders; i += 1)
		folderName = GetIndexedObjNameDFR(dfr, 4, i)
		wList += "root:data:'" + folderName + "':thisWave;"
	endfor
	
	// we need to concatenate these waves into root (all_*)
	String targetWaveList = "Img_VsArea;Img_VsAspectRatio;Img_VsCentre;Img_VsCircularity;Img_VsMajAxis;Img_VsMinAxis;Img_VsPerimeter;Img_MtArea;Img_MtPerimeter;"
	Variable nTargets = ItemsInList(targetWaveList)
	String targetName, tList, conName
	
	SetDataFolder root:
	String fullName,modtList
	
	for(i = 0; i < nTargets; i += 1)
		targetName = StringFromList(i,targetWaveList)
		tList = ReplaceString("thisWave",wList,targetName)
		modtList = tList
		// because some waves might not exist
		for(j = 0; j < numDataFolders; j +=1)
			fullName = StringFromList(j, tList)
			Wave testW = $fullName
			if(!WaveExists(testW))
				modtList = RemoveFromList(fullName,modtList)
			endif
		endfor
		conName = ReplaceString("Img_",targetName,"All_")
		Concatenate/O/NP=0 modtList, $conName
	endfor
	
	// redefine targetList for summation of each folder into row of sum_* wave
	targetWaveList = "Img_VsArea;Img_VsPerimeter;Img_MtArea;Img_MtPerimeter;"
	nTargets = ItemsInList(targetWaveList)
	String wName,tName
	for(i = 0; i < nTargets; i += 1)
		targetName = StringFromList(i,targetWaveList)
		wName = "root:" + ReplaceString("Img_",targetName,"Sum_")
		Make/O/N=(numDataFolders) $wName
		Wave sumW0 = $wName
		tList = ReplaceString("thisWave",wList,targetName)
		for(j = 0; j < numDataFolders; j += 1)
			tName = StringFromList(j,tList)
			Wave tW0 = $tName
			sumW0[j] = sum(tW0)
		endfor
	endfor
	
	// now count the number of Vs and Mt waves in each folder
	Make/O/N=(numDataFolders) Count_Vs,Count_Mt
	targetName = "Img_VsArea"
	tList = ReplaceString("thisWave",wList,targetName)
	for(i = 0; i < numDataFolders; i += 1)
		tName = StringFromList(i,tList)
		Wave tW0 = $tName
		Count_Vs[i] = numpnts(tW0)
	endfor
	targetName = "Img_MtArea"
	tList = ReplaceString("thisWave",wList,targetName)
	for(i = 0; i < numDataFolders; i += 1)
		tName = StringFromList(i,tList)
		Wave tW0 = $tName
		Count_Mt[i] = numpnts(tW0)
	endfor
	WAVE/Z Sum_MtPerimeter,Sum_VsPerimeter
	MatrixOp/O Ratio_CountVsPerCountMito = Count_Vs / Count_Mt
	MatrixOp/O Ratio_CountVsPerSumMitoPerim = Count_Vs / Sum_MtPerimeter
	MatrixOp/O Ratio_SumVsPerimPerSumMitoPerim = Sum_VsPerimeter / Count_Mt
End

// generate the figure
Function MakeTheLayouts()
	DoWindow/K modelLayout
	NewLayout/N=modelLayout
	String modList = WinList("mod_*",";","WIN:1")
	Variable nWindows = ItemsInList(modList)
	Variable PlotsPerPage = 15 // 5 x 3
	String plotName
	String exString = "Tile/A=(" + num2str(ceil(PlotsPerPage/3)) + ",3)"
	
	Variable i,pgNum=1
	
	for(i = 0; i < nWindows; i += 1)
		plotName = StringFromList(i,modList)
		AppendLayoutObject/W=modelLayout/PAGE=(pgnum) graph $plotName
		if(mod((i + 1),PlotsPerPage) == 0 || i == (nWindows -1)) // if page is full or it's the last plot
			LayoutPageAction/W=modelLayout size(-1)=(595, 842), margins(-1)=(18, 18, 18, 18)
			ModifyLayout/W=modelLayout units=0
			ModifyLayout/W=modelLayout frame=0,trans=1
			Execute /Q exString
			if (i != nWindows -1)
				LayoutPageAction/W=modelLayout appendpage
				pgNum += 1
				LayoutPageAction/W=modelLayout page=(pgNum)
			endif
		endif
	endfor
	SavePICT/PGR=(1,-1)/E=-2/W=(0,0,0,0) as "models.pdf"
	
End

//
Function CleanSlate()
	String fullList = WinList("*", ";","WIN:7")
	Variable allItems = ItemsInList(fullList)
	String name
	Variable i
 
	for(i = 0; i < allItems; i += 1)
		name = StringFromList(i, fullList)
		KillWindow/Z $name		
	endfor
	
	KillDataFolder root:data:
		
	// Kill waves in root
//	KillWaves/A/Z
	// Look for data folders and kill them
//	DFREF dfr = GetDataFolderDFR()
//	allItems = CountObjectsDFR(dfr, 4)
//	for(i = 0; i < allItems; i += 1)
//		name = GetIndexedObjNameDFR(dfr, 4, i)
//		KillDataFolder $name		
//	endfor
End