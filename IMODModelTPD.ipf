#pragma TextEncoding = "MacRoman"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

// Menu item for easy execution
Menu "Macros"
	"IMOD Models...",  IMODModelAnalysis()
End

Function IMODModelAnalysis()
	LoadIMODModels()
	SetUpWindows()
	ProcessAllModels()
//	CollectAllMeasurements()
//	MakeSummaryLayout()
//	RotatePits()
//	OverlayAllPits()
	FormatPlots()
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
		// find eigenvectors
		Wave w1 = FindEV(w0)
		TName = CleanupName(currentDF,0) + "_" + NameOfWave(w1)
		AppendToGraph/W=allRotVsPlot w1[][1]/TN=$tName vs w1[][0]
	endfor
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

Function ContourCalcs()
	String wList = WaveList("PM_*",";","")
	String wName = StringFromList(0,wList) // it is probably PM_0
	Wave PMw = $wName
	WAVE/Z distWave, rowPM, pitStartStop
	Variable rStart, rEnd, rFT // rows of PMw
	
	if (pitStartStop[0] < pitStartStop[1])
		rStart = pitStartStop[0]
		rEnd = pitStartStop[1]
	else
		rStart = pitStartStop[1]
		rEnd = pitStartStop[0]
	endif
	
	Make/O/N=1 cPit
	
	Duplicate/O/RMD=[rStart,rEnd][0,2] PMw, pitWave
	cPit[0] = ContourLength(pitWave) // find pit length
	
	Variable nFT = numpnts(rowPM)
	Make/O/N=(nFT) cFT
	
	Variable i
	
	for (i = 0; i < nFT; i += 1)
		rFT = rowPM[i]
		if (rFT < rStart)
			Duplicate/O/RMD=[rFT,rStart][0,2] PMw, FTlocWave
			cFT[i] = ContourLength(FTLocWave) * (-1)
		elseif (rFT == rStart)
			cFT[i] = 0
		else
			Duplicate/O/RMD=[rStart,rFT][0,2] PMw, FTlocWave
			cFT[i] = ContourLength(FTLocWave)
		endif
	endfor
	Duplicate/O cFT, ratioFT
	ratioFT /= cPit[0]
	KillWaves/Z pitWave,FTlocWave,result,result2,result3
End

// Works out the contour length along a line
Function ContourLength(m0)
	Wave m0
	
	Differentiate/METH=1/EP=1/DIM=0 m0
	MatrixOP/O result = m0 * m0
	MatrixOP/O result2 = sumrows(result)
	MatrixOP/O result3 = sqrt(result2)
	Return sum(result3)
End

Function CollectAllMeasurements()
	SetDataFolder root:data:	// relies on earlier load
	DFREF dfr = GetDataFolderDFR()
	String folderName
	Variable numDataFolders = CountObjectsDFR(dfr, 4)
	String wList = ""
	
	Variable i
		
	for(i = 0; i < numDataFolders; i += 1)
		folderName = GetIndexedObjNameDFR(dfr, 4, i)
		wList += "root:data:'" + folderName + "':ratioFT;"
	endfor
	
	SetDataFolder root:
	Concatenate/O/NP=0 wList, allRatioWave
	wList = ReplaceString("ratioFT",wList,"distWave")
	Concatenate/O/NP=0 wList, allDistWave
	Variable nFT = numpnts(allDistWave)
	Make/O/N=(nFT,3) allFTNormWave = 0
	Duplicate/O/FREE allDistWave, w0
	w0 +=50
	w0 /=50	// scaling the distance measurements to a pit that was 100 nm diam
	for(i = 0; i < nFT; i += 1)
		if (allRatioWave[i] >= 0 && allRatioWave[i] <= 1)
			allFTNormWave[i][0] = cos(allRatioWave[i] * PI) * w0[i]
			allFTNormWave[i][1] = sin(allRatioWave[i] * PI) * w0[i]
		elseif (allRatioWave[i] < 0)
			allFTNormWave[i][0] = 1 + (abs(allRatioWave[i]) * PI)
			allFTNormWave[i][1] = w0[i] - 1 // because the distance from membrane to x-axis is 0
		else
			allFTNormWave[i][0] = ((allRatioWave[i] - 1) * -PI) - 1
			allFTNormWave[i][1] = w0[i] - 1 // because the distance from membrane to x-axis is 0
		endif
	endfor
	// now make the mirror version
	Duplicate/O allFTNormWave, allFTNormRWave
	allFTNormRWave[][0] *= -1
	// split them into inside and outside versions
	Duplicate/O allFTNormWave, w1
	Duplicate/O allFTNormRWave, w2
	w1 = (allRatioWave[p] >= 0 && allRatioWave[p] <= 1) ? allFTNormWave[p][q] : NaN
	w2 = (allRatioWave[p] >= 0 && allRatioWave[p] <= 1) ? allFTNormRWave[p][q] : NaN
	Concatenate/O/KILL/NP=0 {w1,w2}, allFTNormIn
	Duplicate/O allFTNormWave, w1
	Duplicate/O allFTNormRWave, w2
	w1 = (allRatioWave[p] >= 0 && allRatioWave[p] <= 1) ? NaN : allFTNormWave[p][q]
	w2 = (allRatioWave[p] >= 0 && allRatioWave[p] <= 1) ? NaN : allFTNormRWave[p][q]
	Concatenate/O/KILL/NP=0 {w1,w2}, allFTNormOut
	// now make waves for scatter plot of distances
	Duplicate/O allDistWave, wIn,wOut
	Concatenate/O/NP=0 {wIn,wOut}, allDistInOut // make a copy of this for allFTNormIn/out plotting
	wIn = (allRatioWave[p] >= 0 && allRatioWave[p] <= 1) ? alldistWave[p] : NaN
	wOut = (allRatioWave[p] >= 0 && allRatioWave[p] <= 1) ? NaN : alldistWave[p]
	WaveTransform zapnans wIn
	WaveTransform zapnans wOut
	Make/O/N=(5,2) distMean = {{-0.1,0.1,NaN,0.9,1.1},{0,0,NaN,0,0}}
	Make/O/N=(5,2) distSD = {{0,0,NaN,1,1},{0,0,NaN,0,0}}
	WaveStats/Q wIn
	distMean[0,1][1] = V_avg
	distSD[0][1] = V_avg - V_sdev
	distSD[1][1] = V_avg + V_sdev
	WaveStats/Q wOut
	distMean[3,4][1] = V_avg
	distSD[3][1] = V_avg - V_sdev
	distSD[4][1] = V_avg + V_sdev
	// need x jitter
	Duplicate/O wIn, xJit
	xJit = 0 + gnoise(0.1)
	Concatenate/O/KILL {xJit,wIn}, distWaveIn
	Duplicate/O wOut, xJit
	xJit = 1 + gnoise(0.1)
	Concatenate/O/KILL {xJit,wOut}, distWaveOut
////	MakeModelCCP()
End


// generate the figure
Function MakeSummaryLayout()
	DoWindow/K summaryLayout
	NewLayout/N=summaryLayout

//	AppendLayoutObject/W=summaryLayout graph pitInPlot
//	AppendLayoutObject/W=summaryLayout graph pitOutPlot
//	AppendLayoutObject/W=summaryLayout graph distPlot

	LayoutPageAction size(-1)=(595, 842), margins(-1)=(18, 18, 18, 18)
	ModifyLayout units=0
	ModifyLayout frame=0,trans=1
//	ModifyLayout left(pitInPlot)=21,top(pitInPlot)=21,width(pitInPlot)=284,height(pitInPlot)=84
//	ModifyLayout left(pitOutPlot)=21,top(pitOutPlot)=110,width(pitOutPlot)=284,height(pitOutPlot)=84
//	ModifyLayout left(distPlot)=428,top(distPlot)=21,width(distPlot)=150,height(distPlot)=150
//	ColorScale/C/N=text0/F=0/A=LT/X=50/Y=1 trace={pitInPlot,allFTNormIn}
//	ColorScale/C/N=text0 "Membrane proximity (nm)"
End
