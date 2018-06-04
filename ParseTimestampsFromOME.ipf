#pragma TextEncoding = "MacRoman"		// For details execute DisplayHelpTopic "The TextEncoding Pragma"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
// These functions will read and parse an OME XML dump from a text file to get timestamps
// Pick Load and parse from Macros menu
// or
// You can do this manually or by using Macros>Load OME XML file
// Then you can pick Parse time stamps from the Macros menu to read time stamps.
// Briefly:
// OME metadata is read in to a single textwave (OMEDumpWave0)
// Then a second function will parse out the image names and numbers,
// creating timestamp waves for each file in an mvd2 library
// Known issue: if the OME-XML file has incomplete timestamps the function makes strange waves

Menu "Macros"
	"Load OME XML file",  LoadTextFileOfOMEXMLdata()
	"Parse time stamps", ParseText()
	"Load and parse", LoadAndParse()
End

Function LoadAndParse()
	LoadTextFileOfOMEXMLdata()
	ParseText()
End

Function LoadTextFileOfOMEXMLdata()
	String pathName	=""	// Name of Igor symbolic path or "" to get dialog
	String fileName	=""	// Name of file to load or "" to get dialog
	
	if ((strlen(pathName)==0) || (strlen(fileName)==0))
		// Display dialog looking for file
		Variable refNum
		String filters = "Text Files (*.txt):.txt;"
		filters += "All Files:.*;"
		Open/D/R/P=$pathName/F=filters/M="Select text file" refNum
		fileName = S_fileName			// S_fileName is set by Open/D
		if (strlen(fileName) == 0)		// User cancelled?
			return -2
		endif
		LoadWave/O/Q/P=$pathName/J/N=OMEDumpWave/K=2 fileName
	endif
End

Function ParseText()
	WAVE/Z/T OMEDumpWave0
	if(!WaveExists(OMEDumpWave0))
		Abort "Please load OME data first."
	endif
	Variable nRows = numpnts(OMEDumpWave0)
	String tval, expr, imageNum, imageName
	Make/O/N=100 imageNumWave,imageRegWave
	Make/O/T/N=100 imageNameWave
	DoWindow/K imageTable
	DoWindow/K matrixTable
	
	Variable i, j = 0
	
	for(i = 0; i < nRows; i += 1)
		tval = OMEDumpWave0[i]
		if(strsearch(tval, "<Image ID=",0) == 0)
			// evaluate string
			expr = "<Image ID=\"Image:([[:digit:]]+)\"."
			SplitString/E=(expr) tval, imageNum
			imageNumWave[j] = str2num(imageNum)
			imageName = ReplaceString("\">",tval,"")	// snip end
			expr = "<Image ID=\"Image:" + imageNum + "\" Name=\""
			imageName = ReplaceString(expr,imageName,"")
			imageNameWave[j] = imageName
			imageRegWave[j] = i
			j += 1
		endif
	endfor
	WaveStats/Q imageNumWave
	DeletePoints V_maxLoc+1, (100-(V_maxLoc+1)), imageNumWave,imageNameWave,imageRegWave
//	Edit /N=imageTable imageNumWave,imageNameWave,imageRegWave

	String deltaT, theC, theT, theZ
	j = 0
	Make/O/N=(nRows,5) mDataWave = NaN
	For(i=0; i<nRows; i+=1)
		tval = OMEDumpWave0[i]
		if(strsearch(tval, "<Plane DeltaT=",0) == 0)
			// evaluate string
			expr = "<Plane DeltaT=\"([[:digit:]\.[:digit:]]+)\"."
			SplitString/E=(expr) tval, deltaT
			mDataWave[j][0] = str2num(deltaT)
			expr = ".TheC=\"([[:digit:]]+)\"."
			SplitString/E=(expr) tval, theC
			mDataWave[j][1] = str2num(theC)
			expr = ".TheT=\"([[:digit:]]+)\"."
			SplitString/E=(expr) tval, theT
			mDataWave[j][2] = str2num(theT)
			expr = ".TheZ=\"([[:digit:]]+)\"."
			SplitString/E=(expr) tval, theZ
			mDataWave[j][3] = str2num(theZ)
			mDataWave[j][4] = i
			j +=1
		EndIf
	EndFor
	WaveStats/Q mDataWave	// max value should be i and should exceed deltaT
	DeletePoints/M=0 V_maxRowLoc+1, (DimSize(mDataWave,0)-V_maxRowLoc), mDataWave
//	Edit/N=matrixTable mDataWave

	// split up 2D wave
	MatrixOp/O c4=col(mDataWave,4)
	nRows = numpnts(c4)	// re-use variable
	Make/O/N=(numpnts(imageNumWave),2) diceWave
	diceWave[0][0] = 0	// set first and last points
	diceWave[numpnts(imageNumWave)-1][1] = c4[nRows-1]
	j = 0
	for(i = 0; i < nRows-1; i += 1)	// i refers to rows in mDataWave
		if((c4[i+1]-c4[i]) > 1)
			diceWave[j][1] = i
			dicewave[j+1][0] = i+1
			j+=1
		endif
	endfor

	// now chop up wave
	String wName
	MatrixOp/O c0 = col(mDataWave,0)
	nRows = dimsize(diceWave,0)
	// how many channels? No longer used
	MatrixOp/O c1 = col(mDataWave,1)
	Variable cVar = wavemax(c1)
	// more than one channel gives same time for each channel 
	Variable multipleVar
	KillWindow/Z timeWaves
	Edit/N=timeWaves
	
	for(i = 0; i < nRows; i += 1)
		if((diceWave[i][1]-diceWave[i][0])>1)
			wName = "dT_" + num2str(imageNumWave[i])
			Duplicate/O/R=[diceWave[i][0],diceWave[i][1]] c0, $wName
			Wave w1 = $wName
			// now get rid of multiples
			multipleVar = WaveChecker(w1)
			if(multipleVar > 0)
				Resample/DOWN=(cVar+1)/N=1 w1
				SetScale/P x 0,1,"", w1
			endif
			AppendToTable/W=timeWaves w1
		endif
	endfor
	KillWaves c4,c1,c0,diceWave
End

STATIC Function WaveChecker(w1)
	Wave w1
	if((sum(w1,pnt2x(w1,0),pnt2x(w1,2))/3) == 0)
		return 3
	elseif((sum(w1,pnt2x(w1,0),pnt2x(w1,1))/2) == 0)
		return 2
	else
		return 0
	endif
End