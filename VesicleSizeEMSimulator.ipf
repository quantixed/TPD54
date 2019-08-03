#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

Function Simulate(rLo,rHi,iter)
	Variable rLo,rHi
	Variable iter
	
	if (rLo >= rHi)
		DoAlert 0, "Check limits"
	endif
	
	Variable step = 5
	Variable simPoints = ceil((rHi - rLo) / step) + 1
	Make/O/N=(simPoints) radiiW = rLo + (p * step)
	Make/O/N=(simPoints) meanW,sdW
	Make/O/N=(iter)/FREE simW
	Variable radius
	
	Variable i,j
	
	for(i = 0; i < simPoints; i += 1)
		radius = radiiW[i]
		for(j = 0; j < iter; j += 1)
			simW[j] = UniformSphere(70, radius)
		endfor
		WaveStats/Q simW
		meanW[i] = V_avg
		sdW[i] = V_sdev
	endfor
	KillWindow/Z p_result
	Display/N=p_result
	AppendToGraph/W=p_result radiiW vs radiiW
	ModifyGraph/W=p_result lstyle(radiiW)=3,rgb(radiiW)=(0,0,0)
	AppendToGraph/W=p_result meanW vs radiiW
	ModifyGraph/W=p_result lsize(meanW)=2
	ErrorBars/W=p_result meanW SHADE= {0,0,(0,0,0,0),(0,0,0,0)},wave=(sdW,sdW)
	Label/W=p_result left "Observed radius (nm)"
	Label/W=p_result bottom "Expected radius (nm)"
	SetAxis/A/N=1/E=1/W=p_result left
	SetAxis/A/N=1/E=1/W=p_result bottom
	ModifyGraph/W=p_result width={Plan,1,bottom,left}
End

///	@param  sectionThickness	thickness of section in nm
///	@param  Radius  radius of sphere in nm
Function UniformSphere(sectionThickness,Radius)
	Variable sectionThickness,Radius
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
	TakeASection(xw,yw,zw,sectionThickness,Radius)
//	Concatenate/O {xw,yw,zw}, triplet
//	NewGizmo
//	AppendToGizmo defaultScatter=triplet
//	ModifyGizmo ModifyObject=scatter0,objectType=scatter,property={ size,0.2}
	if(numpnts(zw) < 100)
		return NaN
	endif
	JointHistogram xw,yw,zw
	WAVE/Z M_JointHistogram
	return 128 * DimDelta(M_JointHistogram,0) + DimOffset(M_JointHistogram,0)
//	// collapse into 2D
//	MatrixOp/O collapseMat = sumbeams(M_JointHistogram)
//	// threshold this view
//	collapseMat[][] = (collapseMat[p][q] > 5) ? collapseMat[p][q] : 0
End

Function TakeASection(xw,yw,zw,thick,Radius)
	Wave xw,yw,zw
	Variable thick // thickness in nm
	Variable Radius // radius of vesicle in nm
	Variable halfThick = thick / 2
	Variable centre = enoise(radius + thick)
	Variable front = centre + (halfThick)
	Variable back = centre - (halfThick)
	zw[] = (zw[p] <= back || zw[p] >= front) ? NaN : zw[p]
	xw[] = (numtype(zw[p]) == 2) ? NaN : xw[p]
	yw[] = (numtype(zw[p]) == 2) ? NaN : yw[p]
	WaveTransform zapnans xw
	WaveTransform zapnans yw
	WaveTransform zapnans zw
End