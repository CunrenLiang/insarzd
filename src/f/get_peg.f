c****************************************************************

      program get_peg

c****************************************************************
c**     
c**   FILE NAME: get_peg.f
c**     
c**   DATE WRITTEN: 6/10/95 
c**     
c**   PROGRAMMER: Scott Hensley
c**     
c**   FUNCTIONAL DESCRIPTION: This program reads simple emphemeris
c**   information and compute the appropriate peg frame as well as 
c**   generating 
c**     
c**   ROUTINES CALLED:none
c**     
c**   NOTES: none
c**     
c**   UPDATE LOG:
c**     This program is modified from get_peg_info.f. Cunren Liang, 21-OCT-2015
c**     Both original program and this program assume WGS84 coordinate
c*****************************************************************

      implicit none


c     PARAMETER STATEMENTS:
      
      real*8 r_awgs84,r_e2wgs84
      parameter(r_awgs84=6378137.d0,r_e2wgs84=.00669437999015d0)

      real*8 pi,r_dtor,r_rtod
      parameter(pi=3.141592653589793d0)     !if you have to ask, give it up
      parameter(r_rtod=180.d0/pi,r_dtor=pi/180.d0)  !radian to degree conversions

      integer i_xyztollh,i_llhtoxyz
      parameter(i_xyztollh=2, i_llhtoxyz=1) 
      integer i_schtoxyz,i_xyztosch
      parameter(i_schtoxyz=0,i_xyztosch=1) 

      character*512 argstring
      integer j,k
      real*8 r_spinvec(3)
      real*8 r_gm, r_earthgm, r_spindot, r_earthspindot


      type ellipsoid
         sequence
         real*8 r_a        
         real*8 r_e2
      end type ellipsoid
      type (ellipsoid) elp

      type peg_struct 
         sequence
         real*8 r_lat
         real*8 r_lon
         real*8 r_hdg
      end type peg_struct
      type (peg_struct) peg

      type pegtrans
         sequence
         real*8 r_mat(3,3)
         real*8 r_matinv(3,3)
         real*8 r_ov(3)
         real*8 r_radcur
      end type pegtrans
      type (pegtrans) ptm

      real*8 r_enumat(3,3),r_xyzenumat(3,3),r_enuvel(3)
      real*8 r_xyzpeg(3),r_vxyzpeg(3),r_llhpeg(3)
      real*8 r_tempv(3), r_tempa(3)
      real*8 r_tempvec(3), r_inertialacc(3), r_bodyacc(3), r_platvel(3)
      real*8 r_platacc(3), r_xyznorm, r_platsch(3)
      real*8 r_xyzschmat(3,3),r_schxyzmat(3,3)
      real*8 r_velnorm


      data r_earthspindot /7.29211573052d-5/
      data r_earthgm /3.98600448073d14/


      elp%r_a = r_awgs84
      elp%r_e2 = r_e2wgs84

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     input the point to be converted: r_xyzpeg, r_vxyzpeg, r_gm, r_spindot
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c      print '(a,$)',' position and velocity of point, Planet GM, Planet Spinrate: x, y, z, vx, vy, vz, r_gm, r_spindot: '
c      read *, r_xyzpeg(1), r_xyzpeg(2), r_xyzpeg(3), r_vxyzpeg(1), r_vxyzpeg(2), r_vxyzpeg(3), r_gm, r_spindot

      if(iargc() .ne. 8)then
         write(6,*) ' '
         write(6,*) ' usage: get_peg x y z vx vy vz gm spin'
         write(6,*) '    x, y, z:    position of platform'
         write(6,*) '    vx, vy, vz: velocity of platform'
         write(6,*) '    gm:         planet gm'
         write(6,*) '    spin:       planet spin'
         write(6,*) ' '
         stop
      endif


      call getarg(1, argstring)
      read(argstring,*) r_xyzpeg(1)

      call getarg(2, argstring)
      read(argstring,*) r_xyzpeg(2)

      call getarg(3, argstring)
      read(argstring,*) r_xyzpeg(3)

      call getarg(4, argstring)
      read(argstring,*) r_vxyzpeg(1)

      call getarg(5, argstring)
      read(argstring,*) r_vxyzpeg(2)

      call getarg(6, argstring)
      read(argstring,*) r_vxyzpeg(3)

      call getarg(7, argstring)
      read(argstring,*) r_gm

      call getarg(8, argstring)
      read(argstring,*) r_spindot


c example:
c     get_peg 3014487.424, 5533954.399, 3400560.970, 3185.529504, 2252.778825, -6468.375282, 398600448073000, 7.29211573052e-05


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     do the calculations now


      write(6,*) ' '
      write(6,'(a,1x,3(f12.3,1x))') 'Pos Peg = ',(r_xyzpeg(j),j=1,3)
      write(6,'(a,1x,3(f12.6,1x))') 'Vel Peg = ',(r_vxyzpeg(j),j=1,3)


      call norm(r_vxyzpeg,r_velnorm)

c     take the lat,lon as the peg point and the heading as the peg heading 

      call latlon(elp,r_xyzpeg,r_llhpeg,i_xyztollh)
      call enubasis(r_llhpeg(1),r_llhpeg(2),r_enumat)
      call tranmat(r_enumat,r_xyzenumat)
      call matvec(r_xyzenumat,r_vxyzpeg,r_enuvel)
      peg%r_hdg = atan2(r_enuvel(1),r_enuvel(2))      

      peg%r_lat = r_llhpeg(1)
      peg%r_lon = r_llhpeg(2) 

      call radar_to_xyz(elp,peg,ptm)

      write(6,*) ' '
      write(6,'(a,1x,f12.7,1x,f12.7,1x,f12.3)') 'Peg Lat/Lon , H = ',
     +     peg%r_lat*r_rtod,peg%r_lon*r_rtod,r_llhpeg(3)
      write(6,'(a,1x,f15.7)') 'Peg Heading = ',peg%r_hdg*r_rtod
      write(6,'(a,1x,f15.5)') 'Radius Curvature = ',ptm%r_radcur

      write(6,*) ' '
      write(6,'(a)') 'Rotation matrix '
      write(6,905) ' First row =  ',ptm%r_mat(1,1),ptm%r_mat(1,2),ptm%r_mat(1,3)
 905  format(a,1x,3(f12.9,1x))
      write(6,905) ' Second row = ',ptm%r_mat(2,1),ptm%r_mat(2,2),ptm%r_mat(2,3)
      write(6,905) ' Third row =  ',ptm%r_mat(3,1),ptm%r_mat(3,2),ptm%r_mat(3,3)
      write(6,*) ' '
      write(6,'(a)') 'Translation vector '
      write(6,906) ' Vector = ',ptm%r_ov
 906  format(a,1x,3(f14.5,1x))

      r_spinvec(1) = 0.
      r_spinvec(2) = 0.
      r_spinvec(3) = r_spindot

      call norm(r_xyzpeg,r_xyznorm)

      call cross(r_spinvec,r_xyzpeg,r_tempv)
      
      do k=1,3
         r_inertialacc(k) = -(r_gm*r_xyzpeg(k))/r_xyznorm**3
      enddo

      call cross(r_spinvec,r_vxyzpeg,r_tempa)
      call cross(r_spinvec,r_tempv,r_tempvec)
      
      do k=1,3
         r_bodyacc(k) = r_inertialacc(k) - 2.d0*r_tempa(k) - r_tempvec(k)
      enddo

c     convert back to a local SCH basis
      
      call convert_sch_to_xyz(ptm,r_platsch,r_xyzpeg,i_xyztosch)
      call schbasis(ptm,r_platsch,r_xyzschmat,r_schxyzmat)
      call matvec(r_xyzschmat,r_bodyacc,r_platacc)
      call matvec(r_xyzschmat,r_vxyzpeg,r_platvel)

      write(6,'(a,x,3(f15.7,x))') 'Platform SCH Velocity (m/s): ',r_platvel
      write(6,'(a,x,3(f15.7,x))') 'Platform SCH Acceleration (m/s^2): ',r_platacc


      end  





