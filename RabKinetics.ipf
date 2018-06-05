#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#include <Waves Average>

Function CoKSAveraging()
	CleanSlate()
	XLLoadWave/S="Sheet1"/R=(A1,PN311)/C=3/W=1/D/K=1 ""
	String wList = WaveList("*",";","")
	String wName
	String expr="([[:alnum:]]+)\\w([[:alpha:]]+)\\w([[:alpha:]]+)\\w([[:digit:]]+)\\w([[:digit:]]+)"
	Variable nWaves = ItemsInList(wList)
	String cond, chan, micy, cell, roi
	Variable cellVar, roiVar
	
	Make/O/N=(nWaves)/T condCellName,condName
	
	Variable i
	
	for(i = 0; i < nWaves; i += 1)
		wName = StringFromList(i,wList)
		SplitString/E=expr wName, cond, chan, micy, cell, roi
		cellVar = str2num(cell)
		roiVar = str2num(roi)
		condCellName[i] = cond + "_" + cell
		condName[i] = cond
	endfor

	// load timestamps
	XLLoadWave/S="Sheet2"/R=(A1,O311)/C=3/W=1/D/K=1 ""
	
	FindDuplicates/RT=uCondCellName condCellname
	FindDuplicates/RT=uCondName condName
	KillWaves condCellName,condName
	
	Variable nCondCell = numpnts(uCondCellName)
	String plotName,srch
	String term = "_G_mito_;_G_cyto_;_R_mito_;_R_cyto_;"
	String tName
	Variable baseline
	Make/O/N=(4,4) colorW = {{0,0,1,1},{1,1,0,0},{0,0,0,0},{0.75,0.25,0.75,0.25}}
	colorW *= 65535

	Variable j
	
	for(i = 0; i < nCondCell; i += 1)
		plotName = "p_" + uCondCellName[i]
		KillWindow/Z $plotName
		Display/HIDE=1/N=$plotName
		tName = "t_" + uCondCellName[i]
		for(j = 0; j < 4; j += 1)
			srch = ReplaceString("_",uCondCellName[i],StringFromList(j,term))
			if(StringMatch(srch,"Rab0_G_mito*") == 1 || StringMatch(srch,"Rab0_G_cyto*") == 1)
				continue
			endif
			wList = WaveList(srch + "_*", ";","")
			Concatenate/O/NP=1/KILL wList, tempMat
			MatrixTranspose TempMat
			MatrixOp/O $srch = averageCols(tempMat)
			KillWaves tempMat
			Wave w0 = $srch
			MatrixTranspose w0
			Redimension/N=-1 w0
			baseline = mean(w0, 0, 5)
			w0 /= baseline
			AppendToGraph/W=$plotName w0 vs $tName
			ModifyGraph/W=$plotName rgb($srch)=(colorW[j][0],colorW[j][1],colorW[j][2],colorW[j][3])
		endfor
		// This is a hack
		// Fake Rab0_G waves
		FakeRab0GWaves()

		TextBox/C/N=text0/F=0/A=LT/X=0.00/Y=0.00 uCondCellName[i]
		ModifyGraph/W=$plotName lsize=2
		Label/W=$plotName left "Fluorescence (F/F\\B0\\M)"
		Label/W=$plotName bottom "Times (s)"
		SetAxis/A/N=1/W=$plotName left
		SetAxis/A/N=1/W=$plotName bottom
	endfor
	MakeTheLayouts()
	DoTheAverages()
End

// generate the figure
Function MakeTheLayouts()
	DoWindow/K summaryLayout
	NewLayout/N=summaryLayout
	String modList = WinList("p_*",";","WIN:1")
	Variable nWindows = ItemsInList(modList)
	Variable PlotsPerPage = 15 // 5 x 3
	String plotName
	String exString = "Tile/A=(" + num2str(ceil(PlotsPerPage/3)) + ",3)"
	
	Variable i,pgNum=1
	
	for(i = 0; i < nWindows; i += 1)
		plotName = StringFromList(i,modList)
		AppendLayoutObject/W=summaryLayout/PAGE=(pgnum) graph $plotName
		if(mod((i + 1),PlotsPerPage) == 0 || i == (nWindows -1)) // if page is full or it's the last plot
			LayoutPageAction/W=summaryLayout size(-1)=(595, 842), margins(-1)=(18, 18, 18, 18)
			ModifyLayout/W=summaryLayout units=0
			ModifyLayout/W=summaryLayout frame=0,trans=1
			Execute /Q exString
			if (i != nWindows -1)
				LayoutPageAction/W=summaryLayout appendpage
				pgNum += 1
				LayoutPageAction/W=summaryLayout page=(pgNum)
			endif
		endif
	endfor
	SavePICT/PGR=(1,-1)/E=-2/W=(0,0,0,0) as "plots.pdf"
End

Function DoTheAverages()
	WAVE/Z/T uCondName
	WAVE/Z colorW
	Variable nRows = numpnts(uCondname)
	String avList,avName,errName
	String term = "_G_mito_;_G_cyto_;_R_mito_;_R_cyto_;"
	String plotName,srch,tList
	
	Variable i,j
	
	for(i = 0; i < nRows; i += 1)
		plotName = "ave_" + uCondName[i]
		KillWindow/Z $plotName
		Display/HIDE=0/N=$plotName
		srch = uCondName[i]
		for(j = 0; j < 4; j += 1)
			srch = uCondName[i] + StringFromList(j,term)
			avlist = Wavelist(srch + "*",";","")
			tList = Wavelist("t_" + uCondName[i] + "*",";","")
			avname = RemoveEnding("W_Ave_" + srch)
			errname = ReplaceString("Ave", avName, "Err")
			fWaveAverage(avList, tList, 3, 1, AvName, ErrName)
			AppendToGraph/W=$plotName $avname
			ErrorBars/W=$plotName $avname SHADE= {0,0,(0,0,0,0),(0,0,0,0)},wave=($ErrName,$ErrName)
			ModifyGraph/W=$plotName rgb($avName)=(colorW[j][0],colorW[j][1],colorW[j][2],colorW[j][3])
			DeletePoints numpnts($avName)-1,1, $AvName,$ErrName
		endfor

		TextBox/C/N=text0/F=0/A=LT/X=0.00/Y=0.00 uCondName[i]
		ModifyGraph/W=$plotName lsize=2
		Label/W=$plotName left "Fluorescence (F/F\\B0\\M)"
		Label/W=$plotName bottom "Times (s)"
		SetAxis/A/N=1/W=$plotName left
		SetAxis/A/N=1/W=$plotName bottom
	endfor
End

Function CleanSlate()
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

Function FakeRab0GWaves()
	String wList = WaveList("t_Rab0*",";","")
	Variable nWaves = ItemsInList(wList)
	if(nWaves ==0)
		return -1
	endif
	String wName, newName
	
	Variable i
	
	for(i = 0; i < nWaves; i += 1)
		wName = StringFromList(i,wList)
		newName = ReplaceString("t_",wName,"")
		newName = ReplaceString("_",newName,"_G_mito_")
		Duplicate/O $wName, $newName
		Wave w1 = $newName
		w1 = 1
		newName = ReplaceString("_G_mito_",newName,"_G_cyto_")
		Duplicate/O $wName, $newName
		Wave w2 = $newName
		w2 = 1
	endfor
End