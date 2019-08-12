#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

Menu "Macros"
	"Flicker Analysis...", /Q, FlickerWorkflow()
End

Function FlickerWorkflow()
	LoadTIFFFilesForAnalysis()
	ShowTheData()
End

Function LoadTIFFFilesForAnalysis()

	String expDiskFolderName, expDataFolderName
	String FileList, ThisFile
	Variable FileLoop, nWaves, i,j,k
	
	NewPath/O/Q/M="Please find disk folder" ExpDiskFolder
	if (V_flag!=0)
		DoAlert 0, "Disk folder error"
		Return -1
	endif
	PathInfo /S ExpDiskFolder
	ExpDiskFolderName=S_path
	FileList=IndexedFile(expDiskFolder,-1,".tif")
	Variable nFiles = ItemsInList(FileList)
	Make/O/N=(nFiles)/T fileName
	Make/O/N=(nFiles) resultW0,resultW1,resultW2,resultW3,resultW4,resultW5
	Variable nRow,nCol,nLayer
	
	NewDataFolder/O/S root:data
	
	for (FileLoop = 0; FileLoop < nFiles; FileLoop += 1)
		ThisFile = StringFromList(FileLoop, FileList)
		fileName[fileLoop] = ThisFile
		expDataFolderName = "img_" + num2str(FileLoop)
		NewDataFolder/O/S $expDataFolderName
		ImageLoad/O/T=tiff/Q/S=0/C=-1/LR3D/P=expDiskFolder/N=lImage ThisFile
		Wave/Z lImage
		resultW0[fileLoop] = ImgVariancePerPixel(lImage)
		resultW1[fileLoop] = ImgVariancePerLayer(lImage)
		// now make a double precision version
		nRow = DimSize(lImage,0)
		nCol = DimSize(lImage,1)
		nLayer = DimSize(lImage,2)
		Make/O/N=(nRow,nCol,nLayer)/D/FREE mImage
		// first normalisation method (divide by mean per frame)
		Wave w0 = MeanPerframe(lImage)
		mImage[][][] = lImage[p][q][r] / w0[r]
		resultW2[fileLoop] = ImgVariancePerPixel(mImage)
		resultW3[fileLoop] = ImgVariancePerLayer(mImage)
		// second normalisation method (divide by max per frame)
		Wave w1 = MaxPerframe(lImage)
		mImage[][][] = lImage[p][q][r] / w1[r]
		resultW4[fileLoop] = ImgVariancePerPixel(mImage)
		resultW5[fileLoop] = ImgVariancePerLayer(mImage)
		KillWaves/Z lImage,w0,w1
		SetDataFolder root:data:
	endfor
	SetDataFolder root:
End

Function ShowTheData()
	Edit/K=0 root:fileName,root:resultW0,root:resultW1,root:resultW2,root:resultW3,root:resultW4,root:resultW5
	WAVE/Z resultW0,resultW1,resultW2,resultW3,resultW4,resultW5
	WAVE/Z/T fileName
	Display/W=(35,45,430,657) resultW2,resultW3 vs fileName
	ModifyGraph tkLblRot(bottom)=90
	SetAxis/A/N=1 left
	ModifyGraph rgb(resultW3)=(0,0,65535)
	Legend/C/N=text0/J/F=0/A=LT/X=0.00/Y=0.00 "\\s(resultW2) per pixel\r\\s(resultW3) per frame"
End

STATIC Function ImgVariancePerPixel(img0)
	Wave img0
	Variable nRow = DimSize(img0,0)
	Variable nCol = DimSize(img0,1)
	Variable nLayer = DimSize(img0,2)
	Make/O/N=(nRow,nCol)/D/FREE result
	
	Variable i,j
	
	for(i = 0; i < nRow; i += 1)
		for(j = 0; j < nCol; j += 1)
			MatrixOp/O/FREE tempW = beam(img0,i,j)
			result[i][j] = Variance(tempW)
		endfor
	endfor
	return mean(result)
End

STATIC Function ImgVariancePerLayer(img0)
	Wave img0
	Variable nRow = DimSize(img0,0)
	Variable nCol = DimSize(img0,1)
	Variable nLayer = DimSize(img0,2)
	Make/O/N=(nLayer)/D/FREE result
	
	Variable i
	
	for(i = 0; i < nLayer; i += 1)
		MatrixOp/O/FREE tempW = layer(img0,i)
		result[i] = Variance(tempW)
	endfor
	return mean(result)
End

STATIC Function/WAVE MeanPerFrame(img0)
	Wave img0
	Variable nRow = DimSize(img0,0)
	Variable nCol = DimSize(img0,1)
	Variable nLayer = DimSize(img0,2)
	Make/O/N=(nLayer)/D/FREE result0
	
	Variable i
	
	for(i = 0; i < nLayer; i += 1)
		MatrixOp/O tempW = layer(img0,i)
		result0[i] = mean(tempW)
	endfor
	return result0
End

STATIC Function/WAVE MaxPerFrame(img0)
	Wave img0
	Variable nRow = DimSize(img0,0)
	Variable nCol = DimSize(img0,1)
	Variable nLayer = DimSize(img0,2)
	Make/O/N=(nLayer)/D/FREE result1
	
	Variable i
	
	for(i = 0; i < nLayer; i += 1)
		MatrixOp/O tempW = layer(img0,i)
		result1[i] = wavemax(tempW)
	endfor
	return result1
End