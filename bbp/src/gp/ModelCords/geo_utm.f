	subroutine geoutm(rlon,rlat,rx,ry,UTM_PROJECTION_ZONE,iway)

c   Convert geodetic longitude and latitude to UTM, and back.
c   Use iway=ILONGLAT2UTM for long/lat to UTM, IUTM2LONGLAT for UTM to lat/long.

        implicit real*8(a-h,o-z)

c   CAMx v2.03
c
c     UTM_GEO performs UTM to geodetic (long/lat) translation, and back.
c
c     This is a Fortran version of the BASIC program "Transverse Mercator
c     Conversion", Copyright 1986, Norman J. Berls (Stefan Musarra, 2/94)
c     Based on algorithm taken from "Map Projections Used by the USGS"
c     by John P. Snyder, Geological Survey Bulletin 1532, USDI.
c
c     Input/Output arguments:
c
c        rlon                  Longitude (deg, negative for West)
c        rlat                  Latitude (deg)
c        rx                    UTM easting (m)
c        ry                    UTM northing (m)
c        UTM_PROJECTION_ZONE   UTM zone
c        iway                  Conversion type
c                              ILONGLAT2UTM=geodetic to UTM
c                              IUTM2LONGLAT=UTM to geodetic

	integer UTM_PROJECTION_ZONE,iway

	logical SUPPRESS_UTM_PROJECTION
	integer ILONGLAT2UTM, IUTM2LONGLAT
	real*8 semimaj, semimin, scfa, north, east

	parameter (SUPPRESS_UTM_PROJECTION = .false.)
	parameter (ILONGLAT2UTM = 0, IUTM2LONGLAT = 1)
	parameter (semimaj=6378206.4d0, semimin=6356583.8d0, scfa=.9996d0)
	parameter (north=0.d0, east=500000.d0)

	PI=4.d0*datan(1.d0)
	TWO_PI=2.d0*PI
	degrad=PI/180.d0
	raddeg=180.d0/PI
	
	if(SUPPRESS_UTM_PROJECTION) then
	  if(iway==ILONGLAT2UTM) then
	    rx=rlon
	    ry=rlat
	  else
	    rlon=rx
	    rlat=ry
	  endif
	  return
	endif

c   Save original parameters

	rlon_save=rlon
	rlat_save=rlat
	rx_save=rx
	ry_save=ry

c   Define parameters of reference ellipsoid

	e2=1.d0-(semimin/semimaj)**2.d0
	e4=e2*e2
	e6=e2*e4
	ep2=e2/(1.d0-e2)

	if(iway==IUTM2LONGLAT) then
	  xx=rx
	  yy=ry
	else
	  dlon=rlon
	  dlat=rlat
	endif

c   Set Zone parameters. cm: central meridian.

	zone=dble(UTM_PROJECTION_ZONE)
	cm=zone*6.d0-183.d0
	cmr=cm*degrad

c   Lat/Lon to UTM conversion

	if(iway==ILONGLAT2UTM) then

	  rlon=degrad*dlon
	  rlat=degrad*dlat

	  delam=dlon-cm
	  if(delam<-180.d0) delam=delam+360.d0
	  if(delam>180.d0) delam=delam-360.d0
	  delam=delam*degrad

	  f1=(1.d0-e2/4.d0-3.d0*e4/64.d0-5.d0*e6/256.d0)*rlat
	  f2=3.d0*e2/8.d0+3.d0*e4/32.d0+45.d0*e6/1024.d0
	  f2=f2*dsin(2.d0*rlat)
	  f3=15.d0*e4/256.d0+45.d0*e6/1024.d0
	  f3=f3*dsin(4.d0*rlat)
	  f4=35.d0*e6/3072.d0
	  f4=f4*dsin(6.d0*rlat)
	  rm=semimaj*(f1-f2+f3-f4)
	  
	  if(dlat==90.d0.or.dlat==-90.d0) then
	    xx=0.d0
	    yy=scfa*rm
	  else
	    rn=semimaj/dsqrt(1.d0-e2*dsin(rlat)**2)
	    t=dtan(rlat)**2
	    c=ep2*dcos(rlat)**2
	    a=dcos(rlat)*delam

	    f1=(1.d0-t+c)*a**3/6.d0
	    f2=5.d0-18.d0*t+t**2+72.d0*c-58.d0*ep2
	    f2=f2*a**5/120.d0
	    xx=scfa*rn*(a+f1+f2)
	    f1=a**2/2.d0
	    f2=5.d0-t+9.d0*c+4.d0*c**2
	    f2=f2*a**4/24.d0
	    f3=61.d0-58.d0*t+t**2+600.d0*c-330.d0*ep2
	    f3=f3*a**6/720.d0
	    yy=scfa*(rm+rn*dtan(rlat)*(f1+f2+f3))
	  endif
	  
	  xx=xx+east
	  yy=yy+north

c   UTM to Lat/Lon conversion

	else

	  xx=xx-east
	  yy=yy-north
	  e1=dsqrt(1.d0-e2)
	  e1=(1.d0-e1)/(1.d0+e1)
	  rm=yy/scfa
	  u=1.d0-e2/4.d0-3.d0*e4/64.d0-5.d0*e6/256.d0
	  u=rm/(semimaj*u)

	  f1=3.d0*e1/2.d0-27.d0*e1**3/32.d0
	  f1=f1*dsin(2.d0*u)
	  f2=21.d0*e1**2/16.d0-55.d0*e1**4/32.d0
	  f2=f2*dsin(4.d0*u)
	  f3=151.d0*e1**3.d0/96.d0
	  f3=f3*dsin(6.d0*u)
	  rlat1=u+f1+f2+f3
	  dlat1=rlat1*raddeg
	  
	  if(dlat1>=90.d0.or.dlat1<=-90.d0) then
	    dlat1=dmin1(dlat1,dble(90.d0))
	    dlat1=dmax1(dlat1,dble(-90.d0))
	    dlon=cm
	  else
	    c1=ep2*dcos(rlat1)**2
	    t1=dtan(rlat1)**2
	    f1=1.d0-e2*dsin(rlat1)**2
	    rn1=semimaj/dsqrt(f1)
	    r1=semimaj*(1.d0-e2)/dsqrt(f1**3)
	    d=xx/(rn1*scfa)

	    f1=rn1*dtan(rlat1)/r1
	    f2=d**2/2.d0
	    f3=5.d0+3.d0*t1+10.d0*c1-4.d0*c1**2-9.d0*ep2
	    f3=f3*d**2*d**2/24.d0
	    f4=61.d0+90.d0*t1+298.d0*c1+45.d0*t1**2-252.d0*ep2-3.d0*c1**2
	    f4=f4*(d**2)**3.d0/720.d0
	    rlat=rlat1-f1*(f2-f3+f4)
	    dlat=rlat*raddeg

	    f1=1.d0+2.d0*t1+c1
	    f1=f1*d**2*d/6.d0
	    f2=5.d0-2.d0*c1+28.d0*t1-3.d0*c1**2+8.d0*ep2+24.d0*t1**2
	    f2=f2*(d**2)**2*d/120.d0
	    rlon=cmr+(d-f1+f2)/dcos(rlat1)
	    dlon=rlon*raddeg
	    if(dlon<-180.d0) dlon=dlon+360.d0
	    if(dlon>180.d0) dlon=dlon-360.d0
	  endif
	  
	endif

	if(iway==IUTM2LONGLAT) then
	  rlon=dlon
	  rlat=dlat
	  rx=rx_save
	  ry=ry_save
	else
	  rx=xx
	  ry=yy
	  rlon=rlon_save
	  rlat=rlat_save
	endif

	return
	end
