      integer maxdata
      parameter (maxdata=3500)
      real*8 eqtime,eqlon
      real*8 data(maxdata,2),time(maxdata),dlon(maxdata),dlat(maxdata)
      integer*4 cycle,pass,select(1),ndata,verbose
      character*256 mission,metafile

      real*4 rlat,rlon,rssh,rtime
      integer*4 cycstart,cycend,passstart,passend
      character*256 comment
      
      nmission=10
      do nm=1,nmission
	read(5,*)comment
	if(comment.eq.'END')stop
        read(5,'(A)')mission
!       write(18,*)mission
	read(5,*)cycstart,cycend
	read(5,*)passstart,passend

      verbose = 2
      select(1) = 0
      call getraw_init(mission,verbose)
!     call getraw_options(16,14)   ! old mss_dtu10 being replaced
      call getraw_options(16,18)   ! new mss_dtu15  jason-2,SARAL, Cryosat2
      call getraw_options(10,02)
      call getraw_options(47,01)
      call getraw_limits(47,-.5d0,.5d0)

      do cycle=cycstart,cycend
              write(18,*)'cycle',cycle
      do pass = passstart,passend
              write(18,*)'pass',pass
        call getraw (cycle,pass,1,maxdata,select,time,dlat,dlon,data,ndata,eqtime,eqlon,metafile)
        write(18,*)'KING',ndata
 	write(18,*)eqtime,eqlon
 	write(18,*)metafile
        if(ndata.gt.1)then
		metafile=metafile(43:53)
 		write(18,*)metafile
        open(16,file=metafile,status='unknown',form='unformatted',access='direct',recl=4*4)
	j=0
        do i = 1,ndata
         if(data(i,1).le.10)then
           rlat=dlat(i)
	   rlon=dlon(i)
	   rssh=data(i,1)
	   rtime=time(i)/(60*60*24)
	   j=j+1
           write(16,rec=j) rlat,rlon,rssh,rtime
         endif
       enddo  !ndata 
       endif
      enddo    !  pass
      enddo    !  cycle
      call getraw_stat(0)
      enddo   ! nmission
      end
