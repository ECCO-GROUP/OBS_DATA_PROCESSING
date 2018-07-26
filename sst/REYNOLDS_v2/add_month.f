      real nine(360*160)
      do ij=1,360*160
      nine(ij)=-9999
      enddo
      open(10,file='SST_monthly_r2_2006',form='unformatted',
     *status='old',
     *access='direct',recl=360*160*4)
      write(10,rec=12)nine
      stop
      end
