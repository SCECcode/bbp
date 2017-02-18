
! --------------------------------------------------------------- TprFrctn
      subroutine tprfrctn (frctnFRONT,frctnBACK,Z,NZ)

c  Apply cosine tapers to first frctnfront*nz points of beginning
*  and last frctnback*nz points of time series array Z.
c  Written by Chuck Mueller, USGS.
c  Modified by D. M. Boore on 8/31/88 to eliminate the use of ZNULL;
c  see FBCTPR_CSM for the original version.

*  Dates:  2/13/90 - if ifront or iback is zero, do not apply a taper.
*         12/12/00 - specify taper in terms of fraction of number of
*                    points rather than percent of number (copied from
*                    fbctpr.for)
!         09/29/12 - Small change to trap for frctn = 0.0

      real Z(*)

      PI = 4.0*ATAN(1.0)
      LZ = NZ*frctnFRONT

      if (lz >= 1) then

        SF = PI/LZ
        do I=1,LZ
          F = 0.5*(1.0-COS(SF*(I-1)))
          Z(I) = Z(I)*F
        end do
        
      end if

      LZ = NZ*frctnback

      if (lz >= 1) then
 
        SF = PI/LZ
        do I=NZ,NZ-LZ+1,-1
          F = 0.5*(1.0-COS(SF*(NZ-I)))
          Z(I) = Z(I)*F
        end do

      end if
      
      return
      end
! --------------------------------------------------------------- TprFrctn

