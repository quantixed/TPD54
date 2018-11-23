#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#include <Waves Average>
#include "ParseTimestampsFromOME"

// Menu item for easy execution
Menu "Macros"
	"Load FRAP Results...",  WorkflowForFRAP()
End

Function WorkflowForFRAP()
	SpecifyConditions()
	NormaliseTheseWaves("frapW*",2,5)
	LoadTimeWaves()	
	RenameAllTheWaves()
	PlotOutNormalisedWavesWithTime()
	FitTheWavesUsingTime()
	MakeTheAverageWavesAndFit()
	SummaryOfFits()
	MakeTheFigure()
End

Function SpecifyConditions()
	CleanSlate()
	String cond0 = "HeLa_GFP_*"
	String cond1 = "HeLa_GFPTPD54_*"
	String cond2 = "HeLa_endoGFPTPD54_*"
	String short0 = "GFP"
	String short1 = "GFPTPD54"
	String short2 = "endoGFPTPD54"
	String hstr = "Specify the names of the FRAP images."
		
	Prompt cond0, "Condition 1"
	Prompt cond1, "Condition 2"
	Prompt cond2, "Condition 3"
	Prompt short0, "Short name 1"
	Prompt short1, "Short name 2"
	Prompt short2, "Short name 3"
	DoPrompt/HELP=hstr "Specify", cond0, cond1, cond2, short0, short1, short2
	
	if (V_flag) 
		return -1
	endif
	Make/O/N=(3)/T condWave={cond0,cond1,cond2}
	Make/O/N=(3)/T shortWave={short0,short1,short2}
	CSVLoader()
End

// This function loads CSVs saved from ImageJ of line profiles
// Makes two two-column waves for each image (2 CSVs for red and green)
Function CSVLoader()
	NewDataFolder/O/S root:data
	
	String expDiskFolderName
	String FileList, ThisFile, newName
	Variable FileLoop, nFiles

	Variable counter=0,i
	 
	NewPath/O/Q/M="Please find disk folder" expDiskFolder
	if (V_flag!=0)
		DoAlert 0, "Disk folder error"
		Return -1
	endif
	PathInfo /S expDiskFolder
	expDiskFolderName = S_path
	FileList = IndexedFile(expDiskFolder,-1,".csv")
	nFiles = ItemsInList(FileList)
	Make/O/N=(nFiles)/T root:csvNameWave
	Wave/T csvNameWave = root:csvNameWave

	for(FileLoop = 0; FileLoop < nFiles; FileLoop += 1)
		ThisFile = StringFromList(FileLoop, FileList)
		// store name of file
//		csvNameWave[fileLoop] = RemoveEnding(ThisFile,".csv")
		csvNameWave[fileLoop] = ThisFile
		LoadWave/Q/A/J/D/W/O/L={0,1,0,0,0}/K=1/P=expDiskFolder ThisFile
		WAVE/Z mean1,mean2,mean3
		if(!WaveExists(mean1) || !WaveExists(mean2) || !WaveExists(mean3))
			DoAlert 0, "Problem with CSV load"
			return -1
		endif
		mean1[] -= mean2[p]
		mean3[] -= mean2[p]
		MatrixOp/O dataWave = mean1 / mean3
		Note/K dataWave, num2str(mean(mean3,2,5))
		newName = "frapW_" + num2str(FileLoop)
		MoveWave dataWave,$("root:"+newName)
		WAVEClear dataWave
		// Kill all imported waves
		KillWaves/A/Z
	endfor
	SetDataFolder root:
End

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
	Variable offset
	
	Variable i
	
	for(i = 0; i < nWaves; i += 1)
		wName = StringFromList(i, wList)
		trimmedName = wName[0,strlen(wName)-3]
		tName = "t_" + trimmedName
		Wave/z w0 = $wName
		WaveStats/Q/M=1 w0
		Wave/z w1 = $tName
		if(!WaveExists(w1))
			Print "Missing", tName
		else
			// set t=0 as lowest point of w0
			offset = w1[V_minRowLoc]
			note/K w1, num2str(offset)
			w1 -= offset
			AppendToGraph/W=p_allNWaves w0 vs w1
			sName = ReplaceString("_n",wName,"_s")
			Wave/z w2 = $sName
			AppendToGraph/W=p_allSWaves w2 vs w1
		endif
	endfor
	Label/W=p_allNWaves left "Fluorescence (Normalized)"
	Label/W=p_allNWaves bottom "Time (s)"
	SetAxis/W=p_allNWaves bottom -4,35
	SetAxis/W=p_allNWaves/A/N=1 left
	Label/W=p_allSWaves left "Fluorescence (Scaled)"
	Label/W=p_allSWaves bottom "Time (s)"
	SetAxis/W=p_allSWaves bottom -4,35
	SetAxis/W=p_allSWaves left 0,1.2
End

STATIC Function MakeTheAverageWavesAndFit()
	WAVE/Z/T condWave, shortWave
	if(!WaveExists(condWave) || !WaveExists(shortWave))
		return -1
	endif
	String yList, xList, avName, errName,srchStr
	srchstr = condWave[0]
	yList = Wavelist(srchStr,";","WIN:p_allSWaves")
	xList = ReplaceString("_s",yList,"")	// note that this will fail with strings containing this substring
	xList = ReplaceString(RemoveEnding(srchStr),xList,"t_" + RemoveEnding(srchStr))
	avName = "W_Ave_" + shortWave[0]
	errName = ReplaceString("Ave", avName, "Err")
	fWaveAverage(yList, xList, 1, 1, AvName, ErrName)
	srchstr = condWave[1]
	yList = Wavelist(srchStr,";","WIN:p_allSWaves")
	xList = ReplaceString("_s",yList,"")	
	xList = ReplaceString(RemoveEnding(srchStr),xList,"t_" + RemoveEnding(srchStr))
	avName = "W_Ave_" + shortWave[1]
	errName = ReplaceString("Ave", avName, "Err")
	fWaveAverage(yList, xList, 1, 1, AvName, ErrName)
	srchstr = condWave[2]
	yList = Wavelist(srchStr,";","WIN:p_allSWaves")
	xList = ReplaceString("_s",yList,"")	
	xList = ReplaceString(RemoveEnding(srchStr),xList,"t_" + RemoveEnding(srchStr))
	avName = "W_Ave_" + shortWave[2]
	errName = ReplaceString("Ave", avName, "Err")
	fWaveAverage(yList, xList, 1, 1, AvName, ErrName)
	// make the plots
	Wave a0 = $("W_Ave_" + shortWave[0])
	Wave a1 = $("W_Ave_" + shortWave[1])
	Wave a2 = $("W_Ave_" + shortWave[2])
	Wave e0 = $("W_Err_" + shortWave[0])
	Wave e1 = $("W_Err_" + shortWave[1])
	Wave e2 = $("W_Err_" + shortWave[2])
	WAVE/Z colorW
	String plotName = "theAveFits"
	Display/N=$plotName a0,a1,a2
	ModifyGraph/W=$plotName lsize=2
	ModifyGraph/W=$plotName rgb($NameOfWave(a0))=(colorW[0][0],colorW[0][1],colorW[0][2])
	ModifyGraph/W=$plotName rgb($NameOfWave(a1))=(colorW[1][0],colorW[1][1],colorW[1][2])
	ModifyGraph/W=$plotName rgb($NameOfWave(a2))=(colorW[2][0],colorW[2][1],colorW[2][2])
	ErrorBars/W=$plotName $NameOfWave(a0) SHADE= {0,0,(0,0,0,0),(0,0,0,0)},wave=(e0,e0)
	ErrorBars/W=$plotName $NameOfWave(a1) SHADE= {0,0,(0,0,0,0),(0,0,0,0)},wave=(e1,e1)
	ErrorBars/W=$plotName $NameOfWave(a2) SHADE= {0,0,(0,0,0,0),(0,0,0,0)},wave=(e2,e2)
	Label/W=$plotName bottom "Time (s)"
	Label/W=$plotName left "Fluorescence (scaled)"
	Legend/W=$plotName/C/N=text0/J/F=0/B=1/A=LT/X=0/Y=0 ""
	Variable lVar,rVar
	lVar = x2pnt(a0,0)
	rVar = numpnts(a0) - 2
	CurveFit dblexp_XOffset a0[lVar,rVar] /W=e0 /I=1 /D 
	lVar = x2pnt(a1,0)
	rVar = numpnts(a1) - 2
	CurveFit dblexp_XOffset a1[lVar,rVar] /W=e1 /I=1 /D 
	lVar = x2pnt(a2,0)
	rVar = numpnts(a2) - 2
	CurveFit dblexp_XOffset a2[lVar,rVar] /W=e2 /I=1 /D
	ModifyGraph/W=$plotName lstyle=3
	ModifyGraph/W=$plotName lstyle($NameOfWave(a0))=0,lstyle($NameOfWave(a1))=0,lstyle($NameOfWave(a2))=0
	SetAxis/W=$plotName bottom -4,35
	SetAxis/W=$plotName left 0,1.2
End

Function LoadTimeWaves()
	SetDataFolder root:
	LoadAndParse()
	WAVE/Z imageRegWave,mDataWave
	WAVE/Z/T OMEDumpWave0
	KillWaves/Z imageRegWave
	KillWaves/Z mDataWave
	KillWaves/Z OMEDumpWave0
End

Function RenameAllTheWaves()
	WAVE/Z/T csvNameWave // list of csv names
	WAVE/Z/T imageNameWave	// list of time waves
	// typically we have csv names like "NameOfMVD2Library - cell 2.csv"
	// time waves have names like "cell 2"
	if(!WaveExists(csvNameWave) || !WaveExists(imageNameWave))
		DoAlert 0, "Problem with the name waves"
		return -1
	endif
	// make integer wave for csvNameWave
	Make/O/N=(numpnts(csvNameWave)) csvNumWave=p
	WAVE/Z imageNumWave	// we already have this one
	if(numpnts(csvNumWave) != numpnts(imageNumWave))
		Print "unequal number of time waves and FRAP waves"
	endif
	Sort csvNameWave,csvNameWave,csvNumWave
	Sort imageNameWave,imageNameWave,imageNumWave
	// at this point we could divine the names automatically
	// to do this we can find the longest substring and remove it
	// then get Igor to find the next longest substrings from sorted lists
	String timeWFileName,timeWName,csvWName,newName,searchString
	Variable timeWNum,csvWNum,theRow
	Variable nRows = numpnts(imageNameWave)
	Variable csvNames = numpnts(csvNumWave)
	Variable i,j
	
	for(i = 0; i < nRows; i += 1)
		timeWFileName = imageNameWave[i]
		timeWNum = imageNumWave[i]
		timeWName = "dT_" + num2str(timeWNum)
		searchString = timeWFileName + ".csv"
		// find substring in csvNameWave
		for(j = 0; j < csvNames; j += 1)
			if(strsearch(csvNameWave[j],searchString,0) != -1)
				theRow = j
				break
			else
				theRow = 999
			endif
		endfor
		if(theRow != 999)
			csvWNum = csvNumWave[theRow]
			csvWName = "frapW_" + num2str(csvWNum)
			// if hyphens are used, remove rather than replace with underscore
			newName = ReplaceString("-",timeWFileName,"")
			// let's call files that end without " n" where n is integer, " 0"
			if(GrepString(newName," [0-9]+") == 0)
				newName += " 0"
			endif
			newName = CleanupName(newName,0)
			// Rename the waves
			Rename $timeWName, $("t_" + newName)
			Rename $csvWName, $newName
			Rename $(csvWName + "_n"), $(newName + "_n")
			Rename $(csvWName + "_s"), $(newName + "_s")
			// for error-checking
		//	Print timeWFileName, "was", timeWName, "renamed as", "t_" + newName, csvNameWave[theRow], "was", csvWName, "renamed as", newName
		else
			Print "No match for", timeWFileName
		endif
	endfor
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
	WAVE/Z/T condWave
	String paramStr = TWave2StringList(condWave)
	ClassifyWaves(fitMat_Name,paramStr)
	LookAtChiSq()
	RecolorTraces("p_allNWaves",0)
	RecolorTraces("p_allSWaves",0)
	RecolorTraces("fit_*",1)
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
	Make/O/N=(3,3) colorW = {{136,221,17},{204,204,119},{238,119,51}}
//	Make/O/N=(3,3) colorW = {{68,136,17},{170,204,119},{153,238,51}}
	colorW *=257
	Make/O/N=(3,4) colorAW
	colorAW[][0,2] = colorW[p][q]
	colorAW[][3] = 32768
End

Function LookAtChiSq()
	WAVE fitMat_chiSq,fitMat_class,colorAW
	KillWindow/Z chiSqPlot
	Display/N=chiSqPlot fitMat_chiSq[][1] vs fitMat_chiSq[][0]
	ModifyGraph/W=chiSqPlot mode=3
	ModifyGraph/W=chiSqPlot zColor(fitMat_chiSq)={fitMat_class,*,*,cindexRGB,0,colorAW}
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
	WAVE/Z/T fitMat_Name,condWave
	WAVE/Z fitMat_Class,colorW,colorAW
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
			if(stringmatch(traceName,condWave[0]) == 1)
				ModifyGraph/W=$graphStr rgb($traceName)=(colorAW[0][0],colorAW[0][1],colorAW[0][2],colorAW[0][3])
			elseif(stringmatch(traceName,condWave[1]) == 1)
				ModifyGraph/W=$graphStr rgb($traceName)=(colorAW[1][0],colorAW[1][1],colorAW[1][2],colorAW[1][3])
			elseif(stringmatch(traceName,condWave[2]) == 1)
				ModifyGraph/W=$graphStr rgb($traceName)=(colorAW[2][0],colorAW[2][1],colorAW[2][2],colorAW[2][3])
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
				if(stringmatch(traceName,condWave[0]) == 1)
					ModifyGraph/W=$plotName rgb($traceName)=(colorW[0][0],colorW[0][1],colorW[0][2])
				elseif(stringmatch(traceName,condWave[1]) == 1)
					ModifyGraph/W=$plotName rgb($traceName)=(colorW[1][0],colorW[1][1],colorW[1][2])
				elseif(stringmatch(traceName,condWave[2]) == 1)
					ModifyGraph/W=$plotName rgb($traceName)=(colorW[2][0],colorW[2][1],colorW[2][2])
				endif
			endfor
		endfor
	endif
End

Function SummaryOfFits()
	WAVE/Z fitMat_double,fitMat_class,colorW,colorAW
	if(!WaveExists(fitmat_double))
		return -1
	endif
	String plotName = "AvsTau"
	KillWindow/Z $plotName
	Display/N=$plotName fitMat_double[][1] vs fitMat_double[][2]
	AppendToGraph/W=$plotName fitMat_double[][3] vs fitMat_double[][4]
	ModifyGraph/W=$plotName mode=3,marker(fitMat_double)=17,marker(fitMat_double#1)=19,zColor(fitMat_double)={fitMat_class,*,*,cindexRGB,0,colorAW},zColor(fitMat_double#1)={fitMat_class,*,*,cindexRGB,0,colorAW}
	SetAxis/W=$plotName left -1,0
	SetAxis/W=$plotName bottom 0,60
	Label/W=$plotName bottom "\\$WMTEX$ \\tau \\$/WMTEX$ (s)"
	Label/W=$plotName left "A"
	Legend/W=$plotName/C/N=text0/J/X=0.00/Y=0.00 "\\s(fitmat_double) Fast\r\\s(fitmat_double#1) Slow"
	// calculate proportions
	Duplicate/O fitmat_double,fitmat_doubleprop
	fitmat_doubleprop[][1] = fitmat_double[p][1] / (fitmat_double[p][1] + fitmat_double[p][3])
	fitmat_doubleprop[][3] = fitmat_double[p][3] / (fitmat_double[p][1] + fitmat_double[p][3])
	plotName = "PropVsTau"
	KillWindow/Z $plotName
	Display/N=$plotName fitMat_doubleprop[][1] vs fitMat_double[][2]
	AppendToGraph/W=$plotName fitMat_doubleprop[][3] vs fitMat_double[][4]
	Label/W=$plotName bottom "\\$WMTEX$ \\tau \\$/WMTEX$ (s)"
	ModifyGraph/W=$plotName mode=3,marker(fitMat_doubleprop)=17,marker(fitMat_doubleprop#1)=19,zColor(fitMat_doubleprop)={fitMat_class,*,*,cindexRGB,0,colorAW},zColor(fitMat_doubleprop#1)={fitMat_class,*,*,cindexRGB,0,colorAW}
	SetAxis/W=$plotName bottom 0,60
	SetAxis/W=$plotName left 0,1
	Label/W=$plotName left "Proportion"
	Legend/W=$plotName/C/N=text0/J/X=0.00/Y=0.00 "\\s(fitmat_doubleprop) Fast\r\\s(fitmat_doubleprop#1) Slow"
	ModifyGraph/W=$plotName mrkThick=0
	// now collect recovery and overall tau into new 2d wave
	WAVE/Z/T fitMat_name
	Variable nWaves = numpnts(fitMat_name)
	Make/O/N=(nWaves,2) fitMat_doubleSummary
	fitMat_doubleSummary[][0] = fitMat_double[p][0]
	String wName, tName
	Variable val50
	Variable i
	for(i = 0; i < nWaves; i += 1)
		wName = fitMat_name[i] + "_dbl"
		Wave w0 = $wName
		Wavestats/Q/M=1 w0
		tName = "t_" + RemoveEnding(fitMat_name[i],"_s")
		Wave t0 = $tName
		// find value where fit crosses 50% of _recovery_
		val50 = ((fitMat_double[i][0] - V_min) / 2) + V_min
		FindLevel/Q w0, val50
		if(V_flag == 1)
			fitMat_doubleSummary[i][1] = NaN
			continue
		endif
		// subtract X0 offset - now that t wave is offset this subtraction is not needed
		fitMat_doubleSummary[i][1] = t0(V_levelX) - t0[V_minRowLoc]
	endfor
	plotName = "RcvryVsOTau"
	KillWindow/Z $plotName
	Display/N=$plotName fitMat_doubleSummary[][0] vs fitMat_doubleSummary[][1]
	Label/W=$plotName bottom "T\B1/2\M (s)"
	Label/W=$plotName left "Recovery"
	ModifyGraph/W=$plotName mode=3,marker=19,zColor(fitMat_doubleSummary)={fitMat_class,*,*,cindexRGB,0,colorAW}
	SetAxis/W=$plotName bottom 0,20
	SetAxis/W=$plotName left 0,1
	ModifyGraph/W=$plotName mrkThick=0
	Wave fitMat_doubleCellF = RetrieveCellFluorescenceFromWaveNotes()
	// make cell intensity plot
	plotName = "IntensityVsOTau"
	KillWindow/Z $plotName
	Display/N=$plotName fitMat_doubleCellF vs fitMat_doubleSummary[][1]
	ModifyGraph/W=$plotName zColor(fitMat_doubleCellF)={fitMat_class,*,*,cindexRGB,0,colorAW}
	ModifyGraph/W=$plotName mode=3,marker=19
	SetAxis/W=$plotName bottom 0,20
	SetAxis/W=$plotName/A/N=1 left
	Label/W=$plotName bottom "T\B1/2\M (s)"
	Label/W=$plotName left "Expression (A.U.)"
	ModifyGraph/W=$plotName mrkThick=0
	// make tau vs tau plot
	plotName = "TauSlowVsFast"
	KillWindow/Z $plotName
	Display/N=$plotName fitMat_double[][4] vs fitMat_double[][2]
	ModifyGraph/W=$plotName mode=3,marker=8
	ModifyGraph/W=$plotName marker=19,zColor(fitMat_double)={fitMat_class,*,*,cindexRGB,0,colorAW}
	SetAxis/W=$plotName left 0,60
	SetAxis/W=$plotName bottom 0,8
	ModifyGraph/W=$plotName zmrkSize(fitMat_double)={fitmat_doubleprop[][3],0.2,1,1,5}
	ModifyGraph/W=$plotName mrkThick=0
	Label/W=$plotName bottom "\\$WMTEX$ \\tau_{fast} \\$/WMTEX$ (s)"
	Label/W=$plotName left "\\$WMTEX$ \\tau_{slow} \\$/WMTEX$ (s)"
	Make/O/N=(3,3) fakeForLegend={{8,8,8},{15,10,5},{0.2,0.5,1.0}}
	AppendToGraph/W=$plotName fakeForLegend[][1] vs fakeForLegend[][0]
	ModifyGraph/W=$plotName mode=3,marker=19,rgb(fakeForLegend)=(34952,34952,34952),zmrkSize(fakeForLegend)={fakeForLegend[*][2],0.2,1,1,5}
	Tag/W=$plotName/C/N=text0/F=0/A=RC/X=-6.00/Y=0.00/L=0 fakeForLegend, 0,"0.2"
	Tag/W=$plotName/C/N=text1/F=0/A=RC/X=-6.00/Y=0.00/L=0 fakeForLegend, 1,"0.5"
	Tag/W=$plotName/C/N=text2/F=0/A=RC/X=-6.00/Y=0.00/L=0 fakeForLegend, 2,"1.0"
End

STATIC Function/WAVE RetrieveCellFluorescenceFromWaveNotes()
	WAVE/Z/T fitMat_name
	Variable nWaves = numpnts(fitMat_Name)
	Make/O/N=(nWaves) fitMat_doubleCellF
	String wName
	
	Variable i
	
	for(i = 0; i < nWaves; i += 1)
		wName = RemoveEnding(fitMat_name[i],"_s")
		Wave w0 = $wName
		fitMat_doubleCellF[i] = str2num(note(w0))
	endfor
	return fitMat_doubleCellF
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

Function MakeTheFigure()
	String layoutName = "finalFigure"
	DoWindow/K $layoutName
	NewLayout/N=$layoutName
	// light touch changes to graphs
	ModifyGraph/W=p_allSWaves noLabel=2
	Label/W=IntensityVsOTau bottom ""
	Legend/W=theAveFits/C/N=text0/J "add legend"
	LayoutPageAction/W=$layoutName size=(595,842),margins=(18,18,18,18)
	AppendLayoutObject/W=$layoutName/F=0/T=1/R=(21,21,258,215) Graph theAveFits
	AppendLayoutObject/W=$layoutName/F=0/T=1/R=(245,21,447,215) Graph RcvryVsOTau
	AppendLayoutObject/W=$layoutName/F=0/T=1/R=(428,21,578,215) Graph TauSlowVsFast
	AppendLayoutObject/W=$layoutName/F=0/T=1/R=(144,115,257,207) Graph p_allSWaves
	AppendLayoutObject/W=$layoutName/F=0/T=1/R=(312,86,447,196) Graph IntensityVsOTau
	String fileName = layoutName + ".pdf"
	SavePICT/O/WIN=$layoutName/PGR=(1,-1)/E=-2/W=(0,0,0,0) as fileName
End

////////////////////////////////////////////////////////////////////////
// Utility functions
////////////////////////////////////////////////////////////////////////
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
	
	// Kill waves in root
	KillWaves/A/Z
	// Look for data folders and kill them
	DFREF dfr = GetDataFolderDFR()
	allItems = CountObjectsDFR(dfr, 4)
	for(i = 0; i < allItems; i += 1)
		name = GetIndexedObjNameDFR(dfr, 4, i)
		if(Stringmatch(name,"*Packages*") != 1)
			KillDataFolder $name		
		endif
	endfor
End

STATIC Function/S TWave2StringList(tW)
	Wave/T tW
	String theString = ""
	Variable nRows = numpnts(tW)
	Variable i
	for(i = 0; i < nRows; i += 1)
		theString += tW[i] + ";"
	endfor
	return theString
End