      SUBROUTINE zib(x,p0,n,p,nx,np0,nn,np,like)

c Zero-inflated Binomial log-likelihood function     

c  Updated 17/01/2007. DH. 

cf2py integer intent(hide),depend(x) :: nx=len(x)
cf2py integer intent(hide),depend(n),check(nn==1 || nn==len(x)) :: nn=len(n)
cf2py integer intent(hide),depend(p),check(np==1 || np==len(x)) :: np=len(p)
cf2py integer intent(hide),depend(p0),check(np0==1 || np0==len(x)) :: np=len(p0)
cf2py double precision intent(out) :: like      
cf2py threadsafe
      IMPLICIT NONE
      INTEGER nx,nn,np,np0,i
      DOUBLE PRECISION like, p(np), p0(np0)
      INTEGER x(nx),n(nn)
      LOGICAL not_scalar_n,not_scalar_p,not_scalar_p0
      INTEGER ntmp
      DOUBLE PRECISION ptmp, p0tmp
      DOUBLE PRECISION factln, thisp
      DOUBLE PRECISION infinity
      PARAMETER (infinity = 1.7976931348623157d308)

      not_scalar_n = (nn .NE. 1)
      not_scalar_p = (np .NE. 1) 
      not_scalar_p0 = (np0 .NE. 1) 

      ntmp = n(1)
      ptmp = p(1)
      p0tmp = p0(1)

      like = 0.0
      do i=1,nx
        if (not_scalar_n) ntmp = n(i)
        if (not_scalar_p) ptmp = p(i)
        if (not_scalar_p0) p0tmp = p0(i)
        
        if ((x(i) .LT. 0) .OR. (ntmp .LT. 0) .OR. (x(i) .GT. ntmp)) then
          like = -infinity
          RETURN
        endif
        
        if ((ptmp .LE. 0.0D0) .OR. (ptmp .GE. 1.0D0)) then
!         if p = 0, number of successes must be 0
          if (ptmp .EQ. 0.0D0) then
            if (x(i) .GT. 0.0D0) then
                like = -infinity
                RETURN
!                 else like = like + 0
            end if
          else if (ptmp .EQ. 1.0D0) then
!           if p = 1, number of successes must be n
            if (x(i) .LT. ntmp) then
                like = -infinity
                RETURN
!                 else like = like + 0
            else if ((x(i).GT.0).OR.(p0tmp.EQ.0.0D0)) then
                like = -infinity
                RETURN
            else
                like = like - dlog(p0tmp)*ntmp
            end if
          else
            like = -infinity
            RETURN
          endif
        else
          thisp = (ptmp**x(i) * (1.-ptmp)**(ntmp-x(i)))
          thisp=thisp*dexp(factln(ntmp)-factln(x(i))-factln(ntmp-x(i)))
          thisp = thisp * (1.-p0tmp)
          if (x(i).EQ.0.0D0) then
            thisp = thisp+p0tmp
          end if
          like = like + dlog(thisp)
        end if
      enddo
      return
      END