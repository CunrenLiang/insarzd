      subroutine intp_coef_norm(nfilter,xintp)

      implicit none
      integer*4     i,j,k,nfilter
      real*4        x,y,pi, sum
      real*4        xintp(0:65544)

      pi = 4.*atan(1.)
c     compute the interpolation factors
      do i=0,nfilter
         j = i*8
         x = real(i)/real(nfilter)
         y = sin(pi*x)/pi
         if(x.ne.0.0  .and. x.ne.1.0) then

            xintp(j  ) = -y/(3.0+x)
            xintp(j+1) =  y/(2.0+x)
            xintp(j+2) = -y/(1.0+x)
            xintp(j+3) =  y/x
            xintp(j+4) =  y/(1.0-x)
            xintp(j+5) = -y/(2.0-x)
            xintp(j+6) =  y/(3.0-x)
            xintp(j+7) = -y/(4.0-x)

c     normalize by the sum of the squares

            sum = 0.
            do k = 0 , 7
               sum = sum + xintp(j+k)**2
            end do
c            sum = sqrt(sum)
            do k = 0 , 7
               xintp(j+k) = xintp(j+k) / sum
            end do

         else if( x.eq.0.0) then
            xintp(j  ) = 0.0
            xintp(j+1) = 0.0
            xintp(j+2) = 0.0
            xintp(j+3) = 1.0
            xintp(j+4) = 0.0
            xintp(j+5) = 0.0
            xintp(j+6) = 0.0
            xintp(j+7) = 0.0
         else if( x.eq.1.0) then
            xintp(j  ) = 0.0
            xintp(j+1) = 0.0
            xintp(j+2) = 0.0
            xintp(j+3) = 0.0
            xintp(j+4) = 1.0
            xintp(j+5) = 0.0
            xintp(j+6) = 0.0
            xintp(j+7) = 0.0
         end if
      end do
      
      return
      end
