#pragma rtGlobals=1		// Use modern global access method.

Function correctedFCS()

	wave decay,macrotimelist,microtimelist,cx

	variable n,tau,dtau,maxt,temp,aa,c,sumdecay,ave,sumb,nt,i
	n=numpnts(macrotimelist)
	maxt=macrotimelist[inf]
	nt=numpnts(cx)
	make/o/n=(nt-1) cor,pwave
	make/o/n=(256,nt-1) bmat=0
	duplicate/o macrotimelist imacrotimelist
	imacrotimelist=maxt-macrotimelist
	reverse imacrotimelist
	duplicate/o microtimelist imicrotimelist
	reverse imicrotimelist
	
	duplicate/o decay a,b
	sumdecay=sum(decay)
	a=(decay>0)? 1 : 0
	ave=sumdecay/sum(a)
	a=decay-a*ave
	a=(decay>0)? a/sqrt(decay) : 0
	aa=norm(a)^2
	
	for(i=0;i<nt-1;i+=1)
		tau=cx[i]
		dtau=cx[i+1]-cx[i]
		calcPhotonAssocDecay(macrotimelist,macrotimelist,microtimelist,n,n,tau,dtau)
		wave W_PhotonAssocDecay
		duplicate/o W_PhotonAssocDecay b1
		sumb=sum(b1)
		calcPhotonAssocDecay(imacrotimelist,imacrotimelist,imicrotimelist,n,n,tau,dtau)
		duplicate/o W_PhotonAssocDecay b2
		b=b1-b2
//		duplicate/o b temp_b
		b=(decay>0)? b/sqrt(decay) : 0
		
//		a=(decay>0)? 1 : 0
//		a=decay-a*ave
//		aa=norm(a)^2
//		c=MatrixDot(a,b)/aa
//		a=(decay>0)? a/sqrt(abs((2*sumb/sumdecay-c)*decay)+abs(c*ave)) : 0
//		aa=norm(a)^2
//		b=(decay>0)? b/sqrt(abs((2*sumb/sumdecay-c)*decay)+abs(c*ave)) : 0
		c=MatrixDot(a,b)/aa
		
		temp=sumb-c*sumdecay			// temp=sumb
		temp/=n*(maxt-tau-dtau/2)/maxt
		temp/=n*dtau/maxt
//		temp_b/=n*(maxt-tau-dtau/2)/maxt
//		temp_b/=n*dtau/maxt
		
		cor[i]=temp
		pwave[i]=c/(maxt-tau-dtau/2)/dtau*maxt
//		bmat[][i]=temp_b[p]
	endfor
	
	killwaves/z imacrotimelist,imicrotimelist
	
end


Function correctedFCS_FLCS()

	wave decay,macrotimeList,microtimeList,cx
	variable i,n,tau,dtau,maxt,temp
	maxt=macrotimelist[inf]
	n=numpnts(macrotimeList)
	make/o/n=(numpnts(cx)-1) cor
	make/o/n=(256,2) em				// make/o/n=(256,1) em
	duplicate/o decay decay_c
	decay_c=(decay>0)? 1 : 0
	em[][1]=1/sum(decay_c)			//
	decay_c=(decay>0)? decay : inf
	temp=wavemin(decay_c)
	decay_c=decay-temp				//
	decay_c=(decay>0)? decay_c : 0
	em[][0]=decay_c/sum(decay_c)
	decay_c=(decay>0)? decay : inf
	matrixop/o filter=Inv(em^T x Inv(Diagonal(decay_c)) x em) x em^T x Inv(Diagonal(decay_c))
	
	for(i=0;i<(numpnts(cx)-1);i+=1)
		tau=cx[i]
		dtau=cx[i+1]-cx[i]
		calc2DCorrHistogram(macrotimeList,macrotimeList,microtimeList,microtimeList,n,n,tau,dtau)
		wave M_CorrHist
		MatrixOP/O filteredM=filter x M_CorrHist x filter^T
		cor[i]=filteredM[0][0]/(dtau*(maxt-tau-dtau/2)/maxt)
	endfor
	
	M_CorrHist=decay[p]*decay[q]
	MatrixOP/O filteredM=filter x M_CorrHist x filter^T
	cor/=filteredM[0][0]/maxt

end
