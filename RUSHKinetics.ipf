#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

// This ipf is not completely automated
// Current method is to do
// 1. NormaliseTheseWaves(searchStr,row0,row1)
// 2. Change t_ waves (time waves from s to min using ConvertSec2Min()
// 3. PlotOutNormalisedWavesWithTime()
// 4. Make the averages - manually
// 5. FitTheWavesUsingTime()

// This function normalises the input waves
// Time waves are prefixes t_* and are excluded
// Searchstr can be "*" for all waves, or "*kd*" etc.
// row0 and row1 specify the start and end of the data to be averaged for baseline
Function NormaliseTheseWaves(searchStr,row0,row1)
	String searchStr
	Variable row0,row1
	
	String wList = WaveList(searchStr,";","")
	String nList = WaveList(searchStr+"_n",";","")
	String sList = WaveList(searchStr+"_s",";","")
	String tList = WaveList("t_*",";","")
	wlist = RemoveFromList(nList,wList)
	wlist = RemoveFromList(sList,wList)
	wList = RemoveFromList(tList,wList)
	
	Variable nWaves = ItemsInList(wList)
	String wName,newRatioName,newScaleName
	Variable bg,mx
	Variable i
	
	// make ratio i.e. F/F0, called *_n - NORMALIZED
	// make scaled wave i.e. 0 to 1, called *_s - SCALED
	for(i = 0; i < nWaves; i += 1)
		wName=StringFromList(i,wList)
		Wave w0 = $wName
		newRatioName = wName + "_n"
		Duplicate/O w0 $newRatioName
		Wave w1 = $newRatioName
		bg = mean(w0,row0,row1)
		w1[] /= bg
		newScaleName = wName + "_s"
		Duplicate /O w0 $newScaleName
		Wave w2 = $newScaleName
		w2[] -= bg
		mx = WaveMax(w2)
		w2[] /= mx
	endfor
End

// this function pairs up the input waves (normalised) with their time counterpart
// note that several input waves share a time wave.
Function PlotOutNormalisedWavesWithTime()
	String wList = WaveList("*_n",";","")
	Variable nWaves = ItemsInList(wList)
	String wName,trimmedName,tName,sName
	KillWindow/Z p_allNWaves
	Display/N=p_allNWaves
	KillWindow/Z p_allSWaves
	Display/N=p_allSWaves
	
	Variable i
	
	for(i = 0; i < nWaves; i += 1)
		wName = StringFromList(i, wList)
		trimmedName = wName[0,strlen(wName)-5]
		if(stringmatch(trimmedName,"X*") == 1)
			trimmedName = ReplaceString("X",trimmedName,"")
		endif
		tName = "t_" + trimmedName
		Wave/z w0 = $wName
		Wave/z w1 = $tName
		if(!WaveExists(w1))
			Print "Missing", tName
		else
			AppendToGraph/W=p_allNWaves w0 vs w1
			sName = ReplaceString("_n",wName,"_s")
			Wave/z w2 = $sName
			AppendToGraph/W=p_allSWaves w2 vs w1
		endif
		if(stringmatch(wName,"*GL*") == 0)
			ModifyGraph/W=p_allNWaves rgb($wName)=(0,0,65535)
			ModifyGraph/W=p_allSWaves rgb($sName)=(0,0,65535)
		endif
	endfor
	Label/W=p_allNWaves left "Ratio (Normalized)"
	Label/W=p_allNWaves bottom "Time (min)"
	Label/W=p_allSWaves left "Ratio (Normalized)"
	Label/W=p_allSWaves bottom "Time (min)"
End

Function FitTheWavesUsingTime()
	String wList = WaveList("*_s",";","")
	// uncomment this line to run on the raw data, also trim 5 characters not 3 below
	wList = ReplaceString("_s",wList,"")
	Variable nWaves = ItemsInList(wList)
	String wName,trimmedName,tName,sName
	Make/O/N=(nWaves)/T fitMat_name
	Make/O/N=(nWaves,4)/D fitMat_a
	Make/O/N=(nWaves,2)/D fitMat_b
	Make/O/N=(nWaves,3)/D resiMat_a=0,resiMat_b=0,halfTMat=0
	Make/O/N=(nWaves)/D fitQW_a=0,fitQW_b=0
	
	Variable i
	
	for(i = 0; i < nWaves; i += 1)
		wName = StringFromList(i, wList)
		trimmedName = wName[0,strlen(wName)-3]
		if(stringmatch(trimmedName,"X*") == 1)
			trimmedName = ReplaceString("X",trimmedName,"")
		endif
		tName = "t_" + trimmedName
		Wave/z w0 = $wName
		Wave/z w1 = $tName
		if(!WaveExists(w1))
			Print "Missing", tName
		else
			fitMat_name[i] = wName
			TheFitter(w0,w1,i)
		endif
	endfor
	SummariseTheFits("GL2","KD")
	MakeTheLayouts("fit_",6,4)
End

// This is a function for fitting
///	@param	yW	the data wave
///	@param	xW	the time wave (X)
///	@param	ii	variable to indicate the iteration that is being called.
Function TheFitter(yW,xW,ii)
	Wave yW,xW
	Variable ii // row number to store W-coefs
	
	WAVE/Z fitMat_a,fitMat_b
	WAVE/Z resiMat_a,resiMat_b
	
	String yName = NameOfWave(yW)
	String plotName = "fit_" + yName
	WaveStats/Q/M=1 yW
	if(V_maxRowLoc == 0)
		return -1
	endif
	KillWindow/Z $plotName
	Display/N=$plotName/HIDE=1 yW vs xW
	if(stringmatch(yName,"*GL*") == 0)
			ModifyGraph/W=$plotName rgb($yName)=(0,0,65535)
	endif
	Variable aRightSide = limit(0,V_maxRowLoc + 2,numpnts(yW)-1)
	String aName = yName + "_a"
	String bName = yName + "_b"
	Make/O/N=(numpnts(yW))/D $aName=NaN,$bName=NaN
	CurveFit/Q hillequation, yW[0,aRightSide] /X=xW /D=$aName
	WAVE/Z W_coef
	fitMat_a[ii][] = W_coef[q]
	AppendToGraph/W=$plotName $aName vs xW
	ModifyGraph/W=$plotName rgb($aName)=(0,0,0)
	// store residuals
	Wave aW = $aName
	ResiStore(yW,aW,resiMat_a,ii)
	if((numpnts(yW) - V_maxRowLoc) > 4)
		CurveFit/Q line, yW[V_maxRowLoc,numpnts(yW)-1] /X=xW /D=$bName
		WAVE/Z W_coef
		fitMat_b[ii][] = W_coef[q]
		AppendToGraph/W=$plotName $bName vs xW
		ModifyGraph/W=$plotName rgb($bName)=(0,0,0)
		Wave bW = $bName
		ResiStore(yW,bW,resiMat_b,ii)
	endif
End

///	@param	yW			This is the data wave
///	@param	fitW		This wave is the fit wave
///	@param	storeW	This is the wave where residuals are stored
///	@param	ii			variable to store iteration
Function ResiStore(yW,fitW,storeW,ii)
	Wave yW,fitW,storeW
	Variable ii
	Duplicate/O/FREE yW,yWcopy
	Duplicate/O/FREE fitW,fitWcopy
	yWcopy[] = (numtype(fitW[p]) == 2) ? NaN : yW[p]
	WaveTransform zapnans yWcopy
	WaveTransform zapnans fitWcopy
	MatrixOp/O/FREE sseW = (yWcopy - fitWcopy) * (yWcopy - fitWcopy)
	storeW[ii][0] = sum(sseW)
	storeW[ii][1] = numpnts(yWcopy)
	storeW[ii][2] = storeW[ii][0] / storeW[ii][1]
End

Function SummariseTheFits(cond0,cond1)
	String cond0,cond1
	
	WAVE/Z/T fitMat_Name
	WAVE/Z fitMat_a,fitMat_b
	WAVE/Z resiMat_a,resiMat_b
	WAVE/Z fitQW_a,fitQW_b // quality waves are 1 if OK, 0 if bad
	// Initial QC - is mean SSE is exceeded?
	fitQW_a[] = (resiMat_a[p][2] > 0.015) ? 0 : 1
	fitQW_b[] = (resiMat_b[p][2] > 0.004) ? 0 : 1
	// Second QC - mark unrealistic values
	fitQW_a[] = (fitMat_a[p][2] > 15 || fitMat_a[p][3] > 40) ? 0 : fitQW_a[p]
	fitQW_b[] = (fitQW_a[p] == 0) ? 0 : fitQW_b[p] // if first fit is bad cancel b
	// Make a temporary version of the fit Matrices and change rows to 0 if fitQW == 0
	Duplicate/O/FREE fitMat_a,tempMat_a
	Duplicate/O/FREE fitMat_b,tempMat_b
	tempMat_a[][] = (fitQW_a[p] == 0) ? 0 : fitMat_a[p][q]
	tempMat_b[][] = (fitQW_b[p] == 0) ? 0 : fitMat_b[p][q]
	// Calculate the half-times
	WAVE/Z halfTMat
	halfTMat[][0] = tempMat_a[p][3] // Thalf for a
	halfTMat[][1] = ((tempMat_a[p][0] + (tempMat_a[p][1] - tempMat_a[p][0]) / 2) - tempMat_b[p][0]) / tempMat_b[p][1]
	halfTMat[][2] = halfTMat[p][1] - halfTMat[p][0] // inter Thalf
	// More QC - any extrapolations longer than 3 h?
	halfTMat[][1] = (halfTMat[p][1] > 150 || halfTMat[p][2] < 6) ? NaN : halfTMat[p][q]
	halfTMat[][2] = (numtype(halfTMat[p][1]) == 2) ? NaN : halfTMat[p][q]
	Variable nRows = numpnts(fitMat_Name)
	Variable count0 = 0
	Variable count1 = 0
	String cellName
	// count instances of each
	Variable i
	for(i = 0; i < nRows; i += 1)
		cellName = fitMat_Name[i]
		if(stringmatch(cellName,"*"+cond0+"*") == 1)
			count0 += 1
		elseif(stringmatch(cellName,"*"+cond1+"*") == 1)
			count1 += 1
		endif
	endfor
	Variable biggestCount = max(count0,count1)
	Make/O/N=(biggestCount,2) fitW_a_50 = NaN, fitW_a_rate = NaN, fitW_b = NaN
	Make/O/N=(biggestCount,2) halfTW_0 = NaN, halfTW_1 = NaN, halfTW_2 = NaN
	// reset counters and store
	count0 = 0
	count1 = 0
	// I don't know why I'm doing it this way -> assignment and zapnans would be easier
	for(i = 0; i < nRows; i += 1)
		cellName = fitMat_Name[i]
		if(stringmatch(cellName,"*"+cond0+"*") == 1)
			fitW_a_50[count0][0] = tempMat_a[i][3] // sigmoid is 2
			fitW_a_rate[count0][0] = tempMat_a[i][2]
			fitW_b[count0][0] = tempMat_b[i][1]
			halfTW_0[count0][0] = halfTMat[i][0]
			halfTW_1[count0][0] = halfTMat[i][1]
			halfTW_2[count0][0] = halfTMat[i][2]
			count0 += 1
		elseif(stringmatch(cellName,"*"+cond1+"*") == 1)
			fitW_a_50[count1][1] = tempMat_a[i][3]
			fitW_a_rate[count1][1] = tempMat_a[i][2]
			fitW_b[count1][1] = tempMat_b[i][1]
			halfTW_0[count1][1] = halfTMat[i][0]
			halfTW_1[count1][1] = halfTMat[i][1]
			halfTW_2[count1][1] = halfTMat[i][2]
			count1 += 1
		endif
	endfor
	// now convert 0 to NaN
	fitW_a_50[][] = (fitW_a_50[p][q] == 0) ? NaN : fitW_a_50[p][q]
	fitW_a_rate[][] = (fitW_a_rate[p][q] == 0) ? NaN : fitW_a_rate[p][q]
	fitW_b[][] = (fitW_b[p][q] == 0) ? NaN : fitW_b[p][q]
	halfTW_0[][] = (halfTW_0[p][q] == 0) ? NaN : halfTW_0[p][q]
	halfTW_1[][] = (halfTW_1[p][q] == 0) ? NaN : halfTW_1[p][q]
	halfTW_2[][] = (halfTW_2[p][q] == 0) ? NaN : halfTW_2[p][q]
	// make labels and then the boxplots
	Make/O/N=2/T labelWave={cond0,cond1}
	MakeBoxPlot("fit_sigmoid50",fitW_a_50,"xHalf (min)")
	MakeBoxPlot("fit_sigmoidRate",fitW_a_rate,"rate")
	MakeBoxPlot("fit_lineRate",fitW_b,"slope (min\S-1\M)")
	MakeBoxPlot("fit_thalf0",halfTW_0,"time (min)")
	MakeBoxPlot("fit_thalf1",halfTW_1,"time (min)")
	MakeBoxPlot("fit_thalf2",halfTW_2,"time (min)")
End

// This function will make a "multicolumn" boxplot (Igor >8 only) 
STATIC Function MakeBoxPlot(plotName,yW,labelString)
	String plotName
	Wave yW
	String labelString
	
	WAVE/Z/T labelWave
	KillWindow/Z $plotName
	Display/N=$plotName
	AppendBoxPlot/W=$plotName yW vs labelWave
	Label/W=$plotName left labelString
	SetAxis/A/N=1/W=$plotName left
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
	String allList = WinList(prefix+"*",";","WIN:1") // edited this line from previous version
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
		plotName = StringFromList(nWindows - 1 - i,modList) // reverse list
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
	String folderStr
	SavePICT/O/WIN=$layoutName/PGR=(1,-1)/E=-2/W=(0,0,0,0) as fileName
End

STATIC Function CleanSlate()
	SetDataFolder root:
	String fullList = WinList("*", ";","WIN:65543")
	Variable allItems = ItemsInList(fullList)
	String name
	Variable i
 
	for(i = 0; i < allItems; i += 1)
		name = StringFromList(i, fullList)
		KillWindow/Z $name		
	endfor
End

Function ConvertSec2Min()
	String wList = WaveList("t_*",";","")
	String wName
	Variable nWaves = ItemsInList(wList)
	Variable i
	for(i = 0; i < nWaves; i += 1)
		wName = StringFromList(i,wList)
		Wave w = $wName
		w /= 60
	endfor
End