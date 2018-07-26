!  Reynolds Temp
!    NOTE:  data need to be filtered
!   output land, missing =0.

      real SST(360,180)
      real ls(360,180)
      character*80 fin,fout

!  land-sea map
      open (9,file='lstags.onedeg.dat',form='unformatted',status='old',
     *access='direct',recl=360*180*4)
      read(9,rec=1)LS
!     write(6,*)ls

!  output file (multiple years, lat limits -79.5 79.5)
      read(5,'(a)')fout
      write(6,'(a)')fout
      open (20,file=fout,form='unformatted',status='unknown',
     *  access='direct',recl=360*160*4)
      IREC=0
      IMO=0

    1 read(5,'(a)')fin
      if(fin.eq.'END')stop
!     write(6,'(a)')fin

      open (10,file=fin,form='unformatted',status='old')

c
c  Loop for each month
  200 READ(10,END=100) 
     *    IYRST,IMST,IDST,IYREND,IMEND,IDEND,NDAYS,INDEX
!     write(6,*)
!    *    IYRST,IMST,IDST,IYREND,IMEND,IDEND,NDAYS,INDEX
      READ (10) SST
      close(10)
      DO  I=1,360
        DO  J=1,180
	  if(LS(i,j).eq.0.)sst(i,j)=0.
      enddo
      enddo

      IMO =IMO + 1
      IREC=IREC+1
      write(20,rec=IREC)((sst(i,j),i=1,360),j=11,170)


c Print date info and SST at one location 
      PRINT 7,IMO,
     * IYRST,IMST,IDST,IYREND,IMEND,IDEND,SST(70+180,80),IREC
    7  FORMAT ('Imonth =',I3,3X,'DATES =',3I4,' - ',3I4,3X,
     * 'SST (110.5W,10.5S) =',F6.2,' IREC =',i5)
       GO TO 1
c Print date info and SST at one location for last week           
  100 PRINT 7,IMO,
     * IYRST,IMST,IDST,IYREND,IMEND,IDEND,SST(70+180,80),IREC

      STOP
      END       

