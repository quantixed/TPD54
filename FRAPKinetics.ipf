#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#include <Waves Average>
// Current method is to do
// 1. NormaliseTheseWaves(searchStr,row0,row1)
// 3. PlotOutNormalisedWavesWithTime()
// 4. Make the averages - manually
// 5. FitTheWavesUsingTime()

// This function normalises the input waves
// Time waves are prefixes t_* and are excluded
// Searchstr can be "*" for all waves, or "*kd*" etc.
// row0 and row1 specify the start and end of the data to be averaged for baseline
// A difference between these procedures and those in RUSH kinetics is that the min and max are reversed
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
	Variable mn,f0
	Variable i
	
	// make ratio i.e. F/F0, called *_n - NORMALIZED
	// make scaled wave i.e. 0 to 1, called *_s - SCALED
	for(i = 0; i < nWaves; i += 1)
		wName=StringFromList(i,wList)
		Wave w0 = $wName
		newRatioName = wName + "_n"
		Duplicate/O w0 $newRatioName
		Wave w1 = $newRatioName
		f0 = mean(w0,row0,row1)
		w1[] /= f0
		newScaleName = wName + "_s"
		Duplicate /O w0 $newScaleName
		Wave w2 = $newScaleName
		mn = WaveMin(w2)
		w2[] -= mn
		f0 = mean(w2,row0,row1)
		w2[] /= f0
	endfor
End

// this function pairs up the input waves (normalised) with their time counterpart
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
		trimmedName = wName[0,strlen(wName)-3]
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
	endfor
	Label/W=p_allNWaves left "Fluorescence (Normalized)"
	Label/W=p_allNWaves bottom "Time (s)"
	Label/W=p_allSWaves left "Fluorescence (Scaled)"
	Label/W=p_allSWaves bottom "Time (s)"
	// Make average waves for scaled
//	String yList, xList, negList, avName, errName
//	yList = Wavelist("GFP_*",";","WIN:p_allSWaves")
//	negList = WaveList("GFP_54*",";","")
//	yList = RemoveFromList(negList,yList)
//	xList = ReplaceString("GFP_",yList,"t_GFP_")
//	avName = "W_Ave_GFP"
//	errName = ReplaceString("Ave", avName, "Err")
//	fWaveAverage(yList, xList, 3, 1, AvName, ErrName)
//	yList = Wavelist("GFP_54*",";","WIN:p_allSWaves")
//	xList = ReplaceString("GFP_",yList,"t_GFP_")
//	avName = "W_Ave_GFP_54"
//	errName = ReplaceString("Ave", avName, "Err")
//	fWaveAverage(yList, xList, 3, 1, AvName, ErrName)
//	yList = Wavelist("endoGFP_*",";","WIN:p_allSWaves")
//	xList = ReplaceString("endoGFP_",yList,"t_endoGFP_")
//	avName = "W_Ave_endoGFP"
//	errName = ReplaceString("Ave", avName, "Err")
//	fWaveAverage(yList, xList, 3, 1, AvName, ErrName)
End

Function FitTheWavesUsingTime()
	String wList = WaveList("*_s",";","")
	Variable nWaves = ItemsInList(wList)
	String wName,trimmedName,tName,sName
	Make/O/N=(nWaves)/T fitMat_name
	Make/O/N=(nWaves,3)/D fitMat_single
	Make/O/N=(nWaves,5)/D fitMat_double
	Make/O/N=(nWaves,2)/D fitMat_chisq
	
	Variable i
	
	for(i = 0; i < nWaves; i += 1)
		wName = StringFromList(i, wList)
		trimmedName = wName[0,strlen(wName)-3]
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
	ClassifyWaves(fitMat_Name,"GFP_*;GFP_54*;endoGFP_*;")
	LookAtChiSq()
	RecolorTraces("p_allNWaves",0)
	RecolorTraces("p_allSWaves",0)
	RecolorTraces("fit_*",1)
//	SummariseTheFits("GL2","KD")
	MakeTheLayouts("fit_",6,4)
End

Function TheFitter(yW,xW,ii)
	Wave yW,xW
	Variable ii // row number to store W-coefs
	
	WAVE/Z fitMat_single,fitMat_double
	WAVE/Z fitMat_chisq
	
	String yName = NameOfWave(yW)
	String plotName = "fit_" + yName
	WaveStats/Q/M=1 yW
	if(V_minRowLoc == 0)
		return -1
	endif
	KillWindow/Z $plotName
	Display/N=$plotName/HIDE=1 yW vs xW
	String sglName = yName + "_sgl"
	String dblName = yName + "_dbl"
	Make/O/N=(numpnts(yW))/D $sglName=NaN,$dblName=NaN
	// single fit to data
	CurveFit/Q exp_XOffset, yW[V_minRowLoc,V_npnts-1] /X=xW /D=$sglName
	// store coef
	WAVE/Z W_coef
	fitMat_single[ii][] = W_coef[q]
	// store chsq
	fitMat_chisq[ii][0] = V_chisq
	AppendToGraph/W=$plotName $sglName vs xW
	ModifyGraph/W=$plotName rgb($sglName)=(32768,0,0)
	
	// double fit to data
	CurveFit/Q dblexp_XOffset, yW[V_minRowLoc,V_npnts-1] /X=xW /D=$dblName
	// store coef
	WAVE/Z W_coef
	fitMat_double[ii][] = W_coef[q]
	// store chsq
	fitMat_chisq[ii][1] = V_chisq
	AppendToGraph/W=$plotName $dblName vs xW
	ModifyGraph/W=$plotName rgb($dblName)=(0,0,0)
End

Function ClassifyWaves(NameWave,condList)
	Wave/T NameWave
	String condList
	Variable nRows = numpnts(NameWave)
	Make/O/N=(nRows) fitMat_class=NaN
	String theString
	Variable i,j
	for(i = 0; i < ItemsInList(condList); i += 1)
		theString = StringFromList(i,condList)
		for(j = 0; j < nRows; j += 1)
			if(stringmatch(NameWave[j],theString) == 1)
				fitMat_class[j] = i
			endif
		endfor
	endfor
	Make/O/N=(3,3) colorW = {{32768,32768,0},{32768,32768,0},{32768,65355,65355}}
End

Function LookAtChiSq()
	WAVE fitMat_chiSq,fitMat_class,colorW
	KillWindow/Z chiSqPlot
	Display/N=chiSqPlot fitMat_chiSq[][1] vs fitMat_chiSq[][0]
	ModifyGraph/W=chiSqPlot mode=3
	ModifyGraph/W=chiSqPlot zColor(fitMat_chiSq)={fitMat_class,*,*,cindexRGB,0,colorW}
	ModifyGraph/W=chiSqPlot width={Aspect,1}
	ModifyGraph/W=chiSqPlot mirror=1;DelayUpdate
	SetAxis/W=chiSqPlot left 0,0.2;DelayUpdate
	SetAxis/W=chiSqPlot bottom 0,0.2
	Label/W=chiSqPlot bottom "Single Exponential (\\$WMTEX$ \\chi^2 \\$/WMTEX$)"
	Label/W=chiSqPlot left "Double Exponential (\\$WMTEX$ \\chi^2 \\$/WMTEX$)"
	SetDrawEnv/W=chiSqPlot dash= 3
	DrawLine/W=chiSqPlot 0,1,1,0
End

Function RecolorTraces(graphStr,optVar)
	String GraphStr
	Variable optVar // 0 for single graph, 1 for many graphs
	WAVE/Z/T fitMat_Name
	WAVE/Z fitMat_Class,colorW
	String traceName, tList, windowList, windowName, plotname
	Variable nTraces, nWindows
	Variable i,j
	
	if(optVar == 0)
		tList = TraceNameList(graphStr,";",1)
		nTraces = ItemsInList(tList)
		for(i = 0; i < nTraces; i += 1)
			traceName = StringFromList(i, tList)
			// this is hardcoded rather than doing lookup
			if(stringmatch(traceName,"*_sgl") == 1 || stringmatch(traceName,"*_dbl") == 1)
				continue
			endif
			if(stringmatch(traceName,"GFP_*") == 1 && stringmatch(traceName,"GFP_54*") == 0)
				ModifyGraph/W=$graphStr rgb($traceName)=(colorW[0][0],colorW[0][1],colorW[0][2])
			elseif(stringmatch(traceName,"GFP_54*") == 1)
				ModifyGraph/W=$graphStr rgb($traceName)=(colorW[1][0],colorW[1][1],colorW[1][2])
			elseif(stringmatch(traceName,"endoGFP_*") == 1)
				ModifyGraph/W=$graphStr rgb($traceName)=(colorW[2][0],colorW[2][1],colorW[2][2])
			endif
		endfor
	elseif(optVar == 1)
		windowList = WinList(graphStr,";","WIN:1")
		nWindows = ItemsInList(windowList)
		for(i = 0; i < nWindows; i += 1)
			plotName = StringFromList(i, windowList)
			tList = TraceNameList(plotName,";",1)
			nTraces = ItemsInList(tList)
			for(j = 0; j < nTraces; j += 1)
				traceName = StringFromList(j, tList)
				// again hard-coded but could be done with lookup
				if(stringmatch(traceName,"*_sgl") == 1 || stringmatch(traceName,"*_dbl") == 1)
					continue
				endif
				if(stringmatch(traceName,"GFP_*") == 1 && stringmatch(traceName,"GFP_54*") == 0)
					ModifyGraph/W=$plotName rgb($traceName)=(colorW[0][0],colorW[0][1],colorW[0][2])
				elseif(stringmatch(traceName,"GFP_54*") == 1)
					ModifyGraph/W=$plotName rgb($traceName)=(colorW[1][0],colorW[1][1],colorW[1][2])
				elseif(stringmatch(traceName,"endoGFP_*") == 1)
					ModifyGraph/W=$plotName rgb($traceName)=(colorW[2][0],colorW[2][1],colorW[2][2])
				endif
			endfor
		endfor
	endif
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

// FRAP movies are <1 min no need to run this
STATIC Function ConvertSec2Min()
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