#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

/// @param	rLo	variable of the smallest vesicle radius in nm that we want to test
/// @param	rHi	variable of the largest vesicle radius in nm that we want to test
/// @param	step	variable of the steps that we will take in nm for the simulation
/// @param	thick	variable for section thickness in nm
/// @param	iter	variable of the number of iterations that we'll do in the simulation
Function Simulate(rLo,rHi,step,thick,iter)
	Variable rLo,rHi,step,thick,iter
	
	if (rLo >= rHi)
		DoAlert 0, "Check limits"
	endif
	
	Variable simPoints = ceil((rHi - rLo) / step) + 1
	Make/O/N=(simPoints) radiiW = rLo + (p * step)
	Make/O/N=(simPoints) meanW,sdW
	Make/O/N=(iter)/FREE simW
	Variable radius
	
	Variable i,j
	
	for(i = 0; i < simPoints; i += 1)
		radius = radiiW[i]
		UniformSphere(radius)
		for(j = 0; j < iter; j += 1)
			simW[j] = ImageTheVesicle(thick, radius)
		endfor
		WaveStats/Q simW
		meanW[i] = V_avg
		sdW[i] = V_sdev
	endfor
	// tidy up
	KillWaves/Z M_JointHistogram
	KillWaves/Z xw,yw,zw
	KillWaves/Z xTemp,yTemp,zTemp
	// make the plot
	PlotTheResults()
End

///	@param  Radius  radius of sphere in nm
STATIC Function UniformSphere(Radius)
	Variable Radius
	
	Variable num = ceil(10 * 4 * pi * Radius^2) // this sets density of points (staining)
	Make/O/N=(num) xw,yw,zw
	Variable phi,theta,rr
	
	Variable i
	
	for(i = 0; i < num; i += 1)
		phi = pi + enoise(pi)
		theta = acos(enoise(1))
		rr = Radius
		xw[i] = rr * sin(theta) * cos(phi)
		yw[i] = rr * sin(theta) * sin(phi)
		zw[i] = rr * cos(theta)
	endfor
	return 1
End

///	@param  sectionThickness	thickness of section in nm
///	@param  Radius  radius of sphere in nm
STATIC Function ImageTheVesicle(sectionThickness,Radius)
	Variable sectionThickness, Radius
	WAVE/Z xw,yw,zw
	Duplicate/O xw,xTemp
	Duplicate/O yw,yTemp
	Duplicate/O zw,zTemp
	TakeASection(xTemp,yTemp,zTemp,sectionThickness,Radius)
	if(numpnts(zTemp) < 100)
		return NaN
	endif
	JointHistogram xTemp,yTemp,zTemp
	WAVE/Z M_JointHistogram
	// return the radius of the object in nm
	return 128 * DimDelta(M_JointHistogram,0) + DimOffset(M_JointHistogram,0)
End

/// @param	xw	wave of x coords
/// @param	yw	wave of y coords
/// @param	zw	wave of z coords
/// @param	thickness	variable of section thickness in nm
/// @param	Radius	variable of vesicle radius in nm
STATIC Function TakeASection(xw,yw,zw,thick,Radius)
	Wave xw,yw,zw
	Variable thick // thickness in nm
	Variable Radius // radius of vesicle in nm
	
	Variable halfThick = thick / 2
	// centre of section may even miss the vesicle
	Variable centre = enoise(radius + thick)
	Variable front = centre + (halfThick)
	Variable back = centre - (halfThick)
	// exclude coords that are outside the section
	zw[] = (zw[p] <= back || zw[p] >= front) ? NaN : zw[p]
	xw[] = (numtype(zw[p]) == 2) ? NaN : xw[p]
	yw[] = (numtype(zw[p]) == 2) ? NaN : yw[p]
	WaveTransform zapnans xw
	WaveTransform zapnans yw
	WaveTransform zapnans zw
End

STATIC Function PlotTheResults()
	String plotName = "p_result"
	WAVE/Z radiiW, meanW, sdW
	Make/O/N=(2) unityW = {0,wavemax(radiiW)}
	KillWindow/Z $plotName
	Display/N=$plotName
	AppendToGraph/W=$plotName unityW vs unityW
	ModifyGraph/W=$plotName lstyle(unityW)=3,rgb(unityW)=(0,0,0)
	AppendToGraph/W=$plotName meanW vs radiiW
	ModifyGraph/W=$plotName lsize(meanW)=2
	ErrorBars/W=$plotName meanW SHADE= {0,0,(0,0,0,0),(0,0,0,0)},wave=(sdW,sdW)
	Label/W=$plotName left "Observed radius (nm)"
	Label/W=$plotName bottom "Expected radius (nm)"
	SetAxis/A/N=1/E=1/W=$plotName left
	SetAxis/A/N=1/E=1/W=$plotName bottom
	ModifyGraph/W=$plotName width={Plan,1,bottom,left}
End