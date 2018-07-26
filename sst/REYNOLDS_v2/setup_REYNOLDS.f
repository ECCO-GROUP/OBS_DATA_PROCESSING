!  Reynolds Temp to Detlef's format
!    NOTE:  data need to be filtered
!   output land, missing =-9999.
!  data in   360 180 -179.5 179.5 -89.5 89.5
!  data out  360 180 .5 359.5 -89.5 89.5

      real SST(360,180)
      real ls(360,180)
      character*80 fin,fout

!  land-sea map
      open (9,file='lstags.onedeg.dat',form='unformatted',status='old',
     *access='direct',recl=360*180*4)
      read(9,rec=1)LS
!     write(6,*)ls

!  output file (multiple years)
      read(5,'(a)')fout
      write(6,'(a)')fout
      open (20,file=fout,form='unformatted',status='unknown',
     *  access='direct',recl=360*180*4)
      IREC=0
      IWK=0

    1 read(5,'(a)')fin
      if(fin.eq.'END')stop
!     write(6,'(a)')fin

      open (10,file=fin,form='unformatted',status='old')

c
c  Loop for each week
  200 READ(10,END=100) 
     *    IYRST,IMST,IDST,IYREND,IMEND,IDEND,NDAYS,INDEX
!     write(6,*)
!    *    IYRST,IMST,IDST,IYREND,IMEND,IDEND,NDAYS,INDEX
      READ (10) SST
      close(10)
      DO  I=1,360
        DO  J=1,180
	  if(LS(i,j).eq.0.)sst(i,j)=-9999.
      enddo
      enddo

      IWK =IWK + 1
      IREC=IREC+1
      write(20,rec=IREC)sst

!     idum=9898
!     i=181
!     xlon=180.5
!     do j = 180,150,-1
!       xlat = float(j) - 90.5
!       print 125, 'lon = ',xlon,'lat = ',xlat,'sst = ',sst(i,j),
!    1     'ice = ',idum,'tagls= ',LS(i,j)
! 125   format(a,f5.1,4x,a,f5.1,4x,a,f5.1,4xa,i3,4x,a,f2.0)
!     enddo

c Print date info and SST at one location 
      PRINT 7,IWK,
     * IYRST,IMST,IDST,IYREND,IMEND,IDEND,SST(70+180,80),IREC
    7  FORMAT ('IWK =',I3,3X,'DATES =',3I4,' - ',3I4,3X,
     * 'SST (110.5W,10.5S) =',F6.2,' IREC =',i5)
       GO TO 1
c Print date info and SST at one location for last week           
  100 PRINT 7,IWK,
     * IYRST,IMST,IDST,IYREND,IMEND,IDEND,SST(70+180,80),IREC

      STOP
      END       

