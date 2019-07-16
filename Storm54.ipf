#pragma TextEncoding = "MacRoman"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

////////////////////////////////////////////////////////////////////////
// Menu items
////////////////////////////////////////////////////////////////////////
Menu "Macros"
	"STORM Spot Analysis...", /Q, DoSTORMAnalysis()
End

////////////////////////////////////////////////////////////////////////
// Master functions and wrappers
////////////////////////////////////////////////////////////////////////

Function DoSTORMAnalysis()
	LoadPrerequisites()
	ClipAndFit()
	CleanOutputs()
	PlotCoefs()
	FindAndPlotWidth()
	MakeTheLayouts("p",5,3)
End

////////////////////////////////////////////////////////////////////////
// Main functions
////////////////////////////////////////////////////////////////////////

Function LoadPrerequisites()
	CleanSlate()
	// Before starting. We need to load in the STORM image
	// Load image of STORM localizations (10 nm/px)
	ImageLoad/O/T=tiff/N=lImage
	// Now we load the coords of maxima from Fiji Find Maxima command
	// Note that Output should have two columns consisting of XY coords in pixels
	LoadWave/A=coord/J/K=1/L={0,1,0,0,0}/O
	WAVE/Z coord0,coord1
	Rename coord0,coordX
	Rename coord1,coordY
End

STATIC Function ClipAndFit()
	Wave/Z coordX = root:coordX
	Wave/Z coordY = root:coordY
	Wave/Z lImage = root:lImage
	if (!WaveExists(lImage) || !WaveExists(coordX) || !WaveExists(coordY))
		DoAlert 0, "Missing wave"
		Return -1
	endif
	
	Variable pxSize = 10 // 10 nm per pixel
	Variable nPoints = numpnts(coordX) // number of maxima
	Variable xMax = dimsize(lImage,0)
	Variable yMax = dimsize(lImage,0)
	Make/O/D/N=(nPoints,7) allCoefs // store the fit coefficients here
	
	Variable xpos, ypos
	String mName
	// for holding cor as zero
	Variable/G K6 = 0
	
	Variable i
	
	for(i = 0; i < nPoints; i += 1)
		xpos = coordX[i]
		ypos = coordY[i]
		// check if the clip will be outside the image and skip it if so
		if(xpos - 20 < 0 || xpos + 20 >= xMax || ypos - 20 < 0 || ypos + 20 >= yMax)
			continue
		endif
		// now make a clip centred on the location 41 x 41 pixels
		mName = "clip_" + num2str(i)
		Duplicate/O/R=[xpos - 20,xpos + 20][ypos - 20,ypos + 20] lImage, $mName
		Wave m0 = $mName
		// do the fit - if debugger finds a singular matrix error, just press go
		CurveFit/Q/H="0000001" Gauss2D m0 /D // hold cor = 0
		WAVE/Z W_coef
		allCoefs[i][] = W_coef[q]
		// the following two lines can be uncommented for debugging or if a stack of TIFFs is needed
//		KillWaves/Z m0
//		KillWaves/Z $("fit_" + mName)
	endfor
	// scale the x and y width of 2D gauss fit to nm
	allCoefs[][3] *= pxsize
	allCoefs[][5] *= pxsize
End

STATIC Function CleanOutputs()
	WAVE/Z allCoefs, coordX, coordY
	if (!WaveExists(allCoefs) || !WaveExists(coordX) || !WaveExists(coordY))
		DoAlert 0, "Missing wave"
		Return -1
	endif
	Variable nPoints = dimsize(allCoefs,0)
	// assemble a "quality wave" and use this to remove spurious fits
	Make/O/FREE/N=(nPoints) qualWave=0
//	// value of 1 means that row has a problem
	// must be between 0-200 nm x or y width
	qualWave[] = (allCoefs[p][3] >= 0 && allCoefs[p][3] <= 200) ? qualWave[p] : 1
	qualWave[] = (allCoefs[p][5] >= 0 && allCoefs[p][5] <= 200) ? qualWave[p] : 1
	// must have a peak that is 0-200 units amplitude
	qualWave[] = (allCoefs[p][1] >= 0 && allCoefs[p][1] <= 200) ? qualWave[p] : 1
	// peak must be within 3 pixels of original location
	qualWave[] = (abs(allCoefs[p][2] - coordX[p]) < 3) ? qualWave[p] : 1
	qualWave[] = (abs(allCoefs[p][4] - coordY[p]) < 3) ? qualWave[p] : 1
	// duplicate then filter
	Duplicate/O allCoefs, allCoefsClean
	allCoefsClean[][] = (qualWave[p] == 0) ? allCoefs[p][q] : NaN
End

STATIC Function PlotCoefs()
	SetDataFolder root:
	WAVE/Z allCoefsClean
	if (!WaveExists(allCoefsClean))
		DoAlert 0, "Missing wave"
		Return -1
	endif
	MatrixOp/O widthX = col(allCoefsClean,3)
	MatrixOp/O widthY = col(allCoefsClean,5)
	MatrixOp/O peakWave = col(allCoefsClean,1)
	// filter out ridiculous values (NaNs inserted by CleanOutputs()
	WaveTransform zapnans widthX
	WaveTransform zapnans widthY
	WaveTransform zapnans peakWave
	Make/O/N=(5,2) fitMean = {{-0.1,0.1,NaN,0.9,1.1},{0,0,NaN,0,0}}
	Make/O/N=(5,2) fitSD = {{0,0,NaN,1,1},{0,0,NaN,0,0}}
	WaveStats/Q widthX
	fitMean[0,1][1] = V_avg
	fitSD[0][1] = V_avg - V_sdev
	fitSD[1][1] = V_avg + V_sdev
	WaveStats/Q widthY
	fitMean[3,4][1] = V_avg
	fitSD[3][1] = V_avg - V_sdev
	fitSD[4][1] = V_avg + V_sdev
	Make/O/N=(2,2) peakMean = {{1.9,2.1},{0,0}}
	Make/O/N=(2,2) peakSD = {{2,2},{0,0}}
	WaveStats/Q peakWave
	peakMean[0,1][1] = V_avg
	peakSD[0][1] = V_avg - V_sdev
	peakSD[1][1] = V_avg + V_sdev
	// need x jitter
	Duplicate/O widthX, xJit
	xJit = 0 + gnoise(0.1)
	Concatenate/O/KILL {xJit,widthX}, fitWaveX
	Duplicate/O widthY, xJit
	xJit = 1 + gnoise(0.1)
	Concatenate/O/KILL {xJit,widthY}, fitWaveY
	Duplicate/O peakWave, xJit
	xJit = 2 + gnoise(0.1)
	Concatenate/O/KILL {xJit,peakWave}, fitWavePeak
	// now plot
	DoWindow/K p_fit
	Display/N=p_fit
	AppendToGraph fitWaveX[][1] vs fitWaveX[][0]
	AppendToGraph fitWaveY[][1] vs fitWaveY[][0]
	ModifyGraph/W=p_fit rgb=(65535,0,0,32768)
	SetAxis/W=p_fit/A/N=1/E=1 left
	Label/W=p_fit left "Fit width (nm)"
	AppendToGraph/R fitWavePeak[][1] vs fitWavePeak[][0]
	ModifyGraph/W=p_fit mode=3,marker=19,msize=2
	ModifyGraph/W=p_fit rgb(fitWavePeak)=(32768,32768,32768,32768)
	ModifyGraph/W=p_fit mrkThick=0
	SetAxis/W=p_fit/A/N=1/E=1 right
	Label/W=p_fit right "Peak (a.u.)"
	SetAxis/W=p_fit bottom -0.5,2.5
	// Hard code the labels
	Make/O/N=3 posWave = p
	Make/O/N=3/T labelWave = {"\u03C3\BX", "\u03C3\BY", "A"}
	ModifyGraph/W=p_fit userticks(bottom)={posWave,labelWave}
	AppendToGraph/W=p_fit fitMean[][1] vs fitMean[][0]
	ModifyGraph/W=p_fit mode(fitMean)=0, lsize(fitMean)=2, rgb(fitMean)=(0,0,0,65535)
	AppendToGraph/W=p_fit fitSD[][1] vs fitSD[][0]
	ModifyGraph/W=p_fit mode(fitSD)=0, lsize(fitSD)=1, rgb(fitSD)=(0,0,0,65535)
	AppendToGraph/W=p_fit/R peakMean[][1] vs peakMean[][0]
	ModifyGraph/W=p_fit mode(peakMean)=0, lsize(peakMean)=2, rgb(peakMean)=(0,0,0,65535)
	AppendToGraph/W=p_fit/R peakSD[][1] vs peakSD[][0]
	ModifyGraph/W=p_fit mode(peakSD)=0, lsize(peakSD)=1, rgb(peakSD)=(0,0,0,65535)
	// Plot xWidth by yWidth
	KillWindow/Z p_fitVs
	Display/N=p_fitVs allCoefsClean[][5] vs allCoefsClean[][3]
	ModifyGraph/W=p_fitVs mode=2
	ModifyGraph/W=p_fitVs log=1
	Label/W=p_fitVs left "\u03C3\BY (nm)"
	Label/W=p_fitVs bottom "\u03C3\BX (nm)"
	SetAxis/W=p_fitVs left 1,200
	SetAxis/W=p_fitVs bottom 1,200
	ModifyGraph/W=p_fitVs width={Aspect,1}
	ModifyGraph/W=p_fitVs grid=1,mirror=1
End

STATIC Function FindAndPlotWidth()
	Wave/Z AllCoefsClean
	Make/O/N=(DimSize(AllCoefsClean,0)) allWidthsClean
	// FWHM is 2.355sigma. Use yWidth only
	allWidthsClean[] = (2 * sqrt(2 * ln(2))) * allCoefsClean[p][5]
	// Histogram from 0 to 100 nm in 1 nm bins
	Make/N=100/O allWidthsClean_Hist
	Histogram/B={0,1,100} allWidthsClean,allWidthsClean_Hist
	KillWindow/Z p_hist
	Display/N=p_hist allWidthsClean_Hist
	ModifyGraph/W=p_hist mode=5,hbFill=0
	Label/W=p_hist bottom "FWHM (nm)"
	SetAxis/W=p_hist bottom 0,100
	SetAxis/W=p_hist/A/N=1 left
	Label/W=p_hist left "Frequency"
	// display stats of the spots in history window
	WaveStats allWidthsClean
	// Make the textbox for the plot
	Duplicate/O/FREE allWidthsClean, tempW
	WaveTransform zapnans tempW
	String meanStr
	sprintf meanStr, "%*.*f nm", 2,2, V_avg
	String medianStr
	sprintf medianStr, "%*.*f nm", 2,2, StatsMedian(tempW)
	String nSpots = num2str(numpnts(tempW))
	TextBox/C/N=text0/F=0/X=0.00/Y=0.00 "Mean = "+meanStr+"\rMedian = "+medianStr+"\rN\Bspots\M = "+nSpots
End

////////////////////////////////////////////////////////////////////////
// Misc functions
////////////////////////////////////////////////////////////////////////

Function MakeStacks()
	// this function will export a TIFF stack of all the clips and all the fits
	// to use, uncomment the killwaves commands as appropriate and run this
	NewPath/Q/O OutputTIFFFolder
	Concatenate/O/KILL/NP=2 WaveList("clip_*",";",""), allTheClips
	Concatenate/O/KILL/NP=2 WaveList("fit_*",";",""), allTheFits
	ImageSave/O/S/U/P=OutputTIFFFolder allTheClips as "allTheClips.tif"
	ImageSave/O/S/U/P=OutputTIFFFolder allTheFits as "allTheFits.tif"
End

////////////////////////////////////////////////////////////////////////
// Utility functions
////////////////////////////////////////////////////////////////////////

STATIC Function CleanSlate()
	String fullList = WinList("*", ";","WIN:7")
	Variable allItems = ItemsInList(fullList)
	String name
	Variable i
 
	for(i = 0; i < allItems; i += 1)
		name = StringFromList(i, fullList)
		KillWindow/Z $name		
	endfor
	
	// Kill waves in root
	KillWaves/A/Z
	// Look for data folders and kill them
	DFREF dfr = GetDataFolderDFR()
	allItems = CountObjectsDFR(dfr, 4)
	for(i = 0; i < allItems; i += 1)
		name = GetIndexedObjNameDFR(dfr, 4, i)
		KillDataFolder $name		
	endfor
End

STATIC Function MakeTheLayouts(prefix,nRow,nCol,[iter, filtVar])
	String prefix
	Variable nRow, nCol
	Variable iter	// this is if we are doing multiple iterations of the same layout
	Variable filtVar // this is the object we want to filter for
	if(ParamIsDefault(filtVar) == 0)
		String filtStr = prefix + "_*_" + num2str(filtVar) + "_*"	// this is if we want to filter for this string from the prefix
	endif
	
	String layoutName = "all"+prefix+"Layout"
	DoWindow/K $layoutName
	NewLayout/N=$layoutName
	String allList = WinList(prefix+"_*",";","WIN:1")
	String modList = allList
	Variable nWindows = ItemsInList(allList)
	String plotName
	
	Variable i
	
	if(ParamIsDefault(filtVar) == 0)
		modList = "" // reinitialise
		for(i = 0; i < nWindows; i += 1)
			plotName = StringFromList(i,allList)
			if(stringmatch(plotName,filtStr) == 1)
				modList += plotName + ";"
			endif
		endfor
	endif
	nWindows = ItemsInList(modList)
	Variable PlotsPerPage = nRow * nCol
	String exString = "Tile/A=(" + num2str(ceil(PlotsPerPage/nCol)) + ","+num2str(nCol)+")"
	
	Variable pgNum=1
	
	for(i = 0; i < nWindows; i += 1)
		plotName = StringFromList(i,modList)
		AppendLayoutObject/W=$layoutName/PAGE=(pgnum) graph $plotName
		if(mod((i + 1),PlotsPerPage) == 0 || i == (nWindows -1)) // if page is full or it's the last plot
			LayoutPageAction/W=$layoutName size(-1)=(595, 842), margins(-1)=(18, 18, 18, 18)
			ModifyLayout/W=$layoutName units=0
			ModifyLayout/W=$layoutName frame=0,trans=1
			Execute /Q exString
			if (i != nWindows -1)
				LayoutPageAction/W=$layoutName appendpage
				pgNum += 1
				LayoutPageAction/W=$layoutName page=(pgNum)
			endif
		endif
	endfor
	String fileName
	if(!ParamIsDefault(iter))
		fileName = layoutName + num2str(iter) + ".pdf"
	else
		fileName = layoutName + ".pdf"
	endif
	if(ParamIsDefault(filtVar) == 0)
		fileName = ReplaceString(".pdf",fileName, "_" + num2str(filtVar) + ".pdf")
	endif
End