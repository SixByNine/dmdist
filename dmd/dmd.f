	program dmd

C  Calls dmdist to get pulsar distance as a function of l, b, and DM.

	real ldeg,l
	real xd(100)
	character*1 limit,xlimit(100)
	common/dxyz/dx0,dy0,dz0
	data rad/57.2957795/

	write(*,1000)
1000	format('At the prompt enter l, b, dm/dist, ndir:'/
     +   ' l:       galactic longitude in deg'/
     +   ' b:       galactic latitude in deg'/
     +   ' dm/dist: dm in cm-3 pc, or dist in kpc'/
     +   ' ndir:    < 0  convert dist to dm'/
     +   '          >= 0 convert dm to dist; ndir > 2 (e.g. ndir=100)',
     +   ' also'/
     +   '               estimates uncertainty in dist.'/)
	dx0=0.
	dy0=0.
	dz0=0.

1	write(*,1001)
1001	format('dmd> ',$)
	read(*,*,err=999,end=999) ldeg,bdeg,dm,ndir
	iters=min(ndir,100)
	if(iters.lt.1) iters=1
	l=ldeg/rad
	b=bdeg/rad

	if(ndir.lt.0) go to 100
	if(ndir.lt.3) go to 20

	idum=0
	do 10 j=1,ITERS
	ds=0.15*gasdev(idum)
	dx0=ds*cos(l)
	dy0=ds*sin(l)
	dz0=0.075*gasdev(idum)
	xdm=dm*(1.0+0.1*gasdev(idum))
	call dmdsm(l,b,1,xdm,xd(j),xlimit(j),sm,smtau,smtheta)
10	continue
	call sort2a(iters,xd,xlimit)
	n1=nint(0.16*iters)
	if(n1.lt.1) n1=1
	n2=iters+1-n1
	dx0=0.
	dy0=0.
	dz0=0.

20	call dmdsm(l,b,1,dm,dist,limit,sm,smtau,smtheta)
	tau=tauiss(dist,smtau,1.0)
	if(iters.lt.3) then
	  write(*,1020) ldeg,bdeg,dm,limit,dist,tau
1020	  format('l:',f6.1,'  b:',f6.1,'  dm:',f7.1,'  dist: ',a1,f5.2,
     +    '  tau: ',f8.3,' ms @ 1GHz')
	else
	  write(*,1030) ldeg,bdeg,xlimit(n1),xd(n1),limit,dist,
     +    xlimit(n2),xd(n2),tau
1030	  format('l:',f6.1,'  b:',f6.1,'  dl: ',a1,f5.2,
     +    '  dist: ',a1,f5.2,'  du: ',a1,f5.2,'  tau: ',f8.3,
     +    ' ms @ 1GHz')
	endif
	go to 1

100	dist=dm
	call dmdsm(l,b,-1,dm,dist,limit,sm,smtau,smtheta)
	tau=tauiss(dist,smtau,1.0)
	write(*,1100) ldeg,bdeg,dist,dm,tau
1100	format('l:',f6.1,'  b:',f6.1,'  dist:',f6.2,
     +    '   dm:',f7.1,'  tau: ',f8.3,' ms @ 1GHz')

	go to 1
	
999	end
