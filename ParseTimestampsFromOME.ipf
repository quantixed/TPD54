#pragma TextEncoding = "MacRoman"		// For details execute DisplayHelpTopic "The TextEncoding Pragma"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
// OME metadata is read in to a single textwave (Textwave0)
// This function will parse out the image names and numbers,
// creating timestamp waves for each file in an mvd2 library
// Save OME-XML output to a file. Load delimited text, pick name Textwave0
// Known issue: if the OME-XML file has incomplete timestamps the function makes strange waves
Function ParseText()
	Wave/T w0 = Textwave0
	Variable nRows = numpnts(w0)
	String tval, expr, imageNum, imageName
	Make/O/N=100 imageNumWave,imageRegWave
	Make/O/T/N=100 imageNameWave
	DoWindow/K imageTable
	DoWindow/K matrixTable
	
	Variable i, j = 0
	
	for(i = 0; i < nRows; i += 1)
		tval = w0[i]
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
	Edit /N=imageTable imageNumWave,imageNameWave,imageRegWave

	String deltaT, theC, theT, theZ
	j = 0
	Make/O/N=(nRows,5) mDataWave = NaN
	For(i=0; i<nRows; i+=1)
		tval = w0[i]
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
	Edit/N=matrixTable mDataWave

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
	// how many channels?
	MatrixOp/O c1 = col(mDataWave,1)
	Variable cVar = wavemax(c1)
	
	for(i = 0; i < nRows; i += 1)
		if((diceWave[i][1]-diceWave[i][0])>1)
			wName = "dT_" + num2str(imageNumWave[i])
			Duplicate/O/R=[diceWave[i][0],diceWave[i][1]] c0, $wName
			Wave w1 = $wName
			// now get rid of multiples
			if((sum(w1,0,cVar)/3)==w1[0])
				Resample/DOWN=(cVar+1)/N=1 w1
			else
				Print "Didn't downsample ", wName // report this case
			endif
		endif
	endfor
	KillWaves c4,c1,c0
End