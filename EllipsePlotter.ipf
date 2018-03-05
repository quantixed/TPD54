#pragma TextEncoding = "MacRoman"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

// suffixList is a string containing "_Ctrl" etc. in semi-colon separated list
Function EllipsePlotter(suffixList)
	String suffixList
	Variable nGroups = ItemsInList(suffixList)
	Make/O/N=81/FREE theta,ellipseY,ellipseX
	theta = 2 * pi * p/80
	Variable xMed,xQ25,xQ75,yMed,yQ25,yQ75
	String wName0, wName1, wName, plotName
	KillWindow/Z summaryLayout
	NewLayout/N=summaryLayout
	LayoutPageAction/W=summaryLayout size(-1)=(595, 842), margins(-1)=(18, 18, 18, 18)
	ModifyLayout units=0
	String exString = "Tile/A=(7,5) "
	
	Variable i
	
	for(i = 0; i < nGroups; i += 1)
		wName0 = "All_VsMajAxis_" + StringfromList(i,suffixList)
		wName1 = "All_VsMinAxis_" + StringfromList(i,suffixList)
		Wave w0 = $wName0
		Wave w1 = $wName1
		StatsQuantiles/Q w0
		xMed = V_Median
		xQ25 = V_Q25
		xQ75 = V_Q75
		StatsQuantiles/Q w1
		yMed = V_Median
		yQ25 = V_Q25
		yQ75 = V_Q75
		plotName = "stats_" + StringfromList(i,suffixList)
		KillWindow/Z $plotName
		Display/N=$plotName
		// Median Wave
		wName = "Vs_Median_" + StringfromList(i,suffixList)
		ellipseY = yMed * cos(theta)
		ellipseX = xMed * sin(theta)
		Concatenate/O {ellipseX,ellipseY}, $wName
		AppendToGraph/W=$plotName $wName[][1] vs $wName[][0]
		ModifyGraph/W=$plotName rgb($wName)=(65535,0,65535)
		// Q25 Wave
		wName = "Vs_Q25_" + StringfromList(i,suffixList)
		ellipseY = yQ25 * cos(theta)
		ellipseX = xQ25 * sin(theta)
		Concatenate/O {ellipseX,ellipseY}, $wName
		AppendToGraph/W=$plotName $wName[][1] vs $wName[][0]
		ModifyGraph/W=$plotName rgb($wName)=(65535,0,65535,16383)
		// Q75 Wave
		wName = "Vs_Q75_" + StringfromList(i,suffixList)
		ellipseY = yQ75 * cos(theta)
		ellipseX = xQ75 * sin(theta)
		Concatenate/O {ellipseX,ellipseY}, $wName
		AppendToGraph/W=$plotName $wName[][1] vs $wName[][0]
		ModifyGraph/W=$plotName rgb($wName)=(65535,0,65535,16383)
		ModifyGraph/W=$plotName lsize=2
		SetAxis/W=$plotName left -50,50
		SetAxis/W=$plotName bottom -50,50
		ModifyGraph/W=$plotName width={Plan,1,bottom,left}
		ModifyGraph/W=$plotName mirror=1
		ModifyGraph/W=$plotName zero=4,standoff=0
		ModifyGraph/W=$plotName margin=26
		AppendLayoutObject/W=summaryLayout graph $plotName
		exString += plotName + ","
	endfor
	ModifyLayout/W=summaryLayout frame=0,trans=1
	DoWindow/F summaryLayout
	Execute/Q exString
End