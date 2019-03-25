c     *****************************************************************
	subroutine aprod(transa, m, n, x, y, rw, iw)
c     Multiply a matrix by a vector
c     Version for use with sparse matrix specified by
c     output of subroutine sparse for use with LSQR

c     original routine from hypoDD by Felix Waldhauser
c     modified by Haijiang Zhang
c     modified by C. Muffels, SSPA, for use with PROPACK and PEST

	implicit none

      integer          :: m,n,i,j,ierr
      integer          :: iw(*)
      double precision :: rw(*),x(*),y(*)
      character*1      :: transa
      
	! iw[1]  Number of non-zero elements in a
	! iw[2:iw[1]+1]  Row indices of non-zero elements
	! iw[iw[1]+2:2*iw[1]+1]  Column indices
      
      do i=1,iw(1)
	  if(transa.eq.'N'.or.transa.eq.'n') then
c	    compute  y = y + a*x
	    y(iw(1+i))=y(iw(1+i))+rw(i)*x(iw(iw(1)+1+i))
      else
c	    compute  y = y + a(transpose)*x
	    y(iw(iw(1)+1+i))=y(iw(iw(1)+1+i))+rw(i)*x(iw(1+i))
	  endif
      end do
      
      return
        
	end subroutine aprod


c     *****************************************************************
      subroutine dzero(n, x , incx)
      implicit none
      integer n, incx
      double precision x(*),r8zero
      parameter (r8zero = 0.0d0)
      integer i

      if ((n.gt.0).and.(incx.ne.0)) then
         if (incx.eq.1) then
            do i=1,n
               x(i) = r8zero
            enddo
         else
            do i=1,n
               x(1+(i-1)*incx) = r8zero
            enddo
         endif
      endif
      return
      end

c     *****************************************************************
      subroutine izero(n, x , incx)
      implicit none
      integer n, incx
      integer x(*),i4zero
      parameter (i4zero = 0)
      integer i

      if ((n.gt.0).and.(incx.ne.0)) then
         if (incx.eq.1) then
            do i=1,n
               x(i) = i4zero
            enddo
         else
            do i=1,n
               x(1+(i-1)*incx) = i4zero
            enddo
         endif
      endif
      return
      end

c     *****************************************************************
      double precision function pdnrm2(n, x, incx)
      implicit none
      integer n, incx
      double precision x(*)
      
      integer i
      double precision sum

      if ((n.gt.0).and.(incx.ne.0)) then    
         sum = 0d0
         if (incx.eq.1) then
c$OMP PARALLEL DO  reduction(+:sum) schedule(static)
            do i=1,n
               sum = sum + x(i)**2
            enddo
         else
c$OMP PARALLEL DO firstprivate(incx) reduction(+:sum) schedule(static)
            do i=1,n
               sum = sum + x(1+(i-1)*incx)**2
            enddo
         endif
         pdnrm2 = sqrt(sum)
      else
         pdnrm2 = 0d0
      endif   
      return
      end

c     *****************************************************************
      subroutine dbdqr(ignorelast, jobq, n, D, E, c1, c2, Qt, ldq)
      implicit none

c Compute QR factorization B = Q*R of (n+1) x n lower bidiagonal matrix 
c with diagonal elements d(1)...d(n) and first subdiagonal elements
c e(1)...e(n). On return [0 ... 0 c1 c2]' = Q'*[0 ... 0 1]'.
c If ignorelast.eq..true. then e(n) is assumed to be zero.
c
c If jobq=='Y' then on return Qt contains Q^T.

c     %------------%
c     | Parameters |
c     %------------%
      character*1 jobq
      logical ignorelast
      integer n,ldq
      double precision D(*),E(*),c1,c2,Qt(ldq,*)
      
c     %-----------------%
c     | Local variables |
c     %-----------------%
      integer i,j
      double precision cs,sn,r

c     %------------------------------------%
c     | External Functions and Subroutines |
c     %------------------------------------%
      logical lsame
      external lsame

c-------------------- Here begins executable code ---------------------

      if (n.lt.1) return
      if (lsame(jobq,'Y')) then
         do j=1,n+1
            do i=1,n+1
               Qt(i,j) = 0.0
            enddo
            Qt(j,j) = 1.0
         enddo
      endif
      do i=1,n-1
         call dlartg(d(i),e(i),cs,sn,r)
         d(i) = r
         e(i) = sn*d(i+1)
         d(i+1) = cs*d(i+1)
         if (lsame(jobq,'Y')) then
            do j=1,i
               Qt(i+1,j) = -sn*Qt(i,j)
               Qt(i,j) = cs*Qt(i,j)
            enddo
            Qt(i,i+1) = sn
            Qt(i+1,i+1) = cs
         endif
      enddo
      if (.not.ignorelast) then
         call dlartg(d(n),e(n),cs,sn,r)
         d(n) = r
         e(n) = 0.0
         c1 = sn
         c2 = cs
         if (lsame(jobq,'Y')) then
            do j=1,i
               Qt(i+1,j) = -sn*Qt(i,j)
               Qt(i,j) = cs*Qt(i,j)
            enddo
            Qt(i,i+1) = sn
            Qt(i+1,i+1) = cs
         endif
      endif
      end

c     *****************************************************************
      subroutine drefinebounds(n,k,theta,bound,tol,eps34)
c
c     Refine Lanczos error bounds using the gap theorem.
c     
c     Input arguments: 
c              n:     smallest dimension of original matrix
c              k:     number of Ritz values to refine
c              theta: array of Ritz values
c              bound: array of unrefined error bounds
c              tol:   clustering tolerance
c              eps34: machine epsilon to the power 3/4.

c     %-----------%
c     | Arguments |
c     %-----------%
      implicit none
      integer n,k
      double precision theta(*), bound(*), tol, eps34

c     %-----------------%
c     | Local variables |
c     %-----------------%
      double precision gap
      integer i,l

c     %------------------------------------%
c     | External Functions and Subroutines |
c     %------------------------------------%
      double precision dlapy2
      external dlapy2

c-------------------- Here begins executable code ---------------------
      if (k.le.1) return
      do i=1,k
         do l=-1,1,2
            if ((l.eq.1.and.i.lt.k) .or. (l.eq.-1.and.i.gt.1)) then
               if (abs(theta(i)-theta(i+l)) .lt. eps34*(theta(i))) then
                  if (bound(i).gt.tol .and. bound(i+l).gt.tol) then
                     bound(i+l) = dlapy2(bound(i),bound(i+l))
                     bound(i) = 0.0
                  endif
               endif
            endif
         enddo         
      enddo
      do i=1,k
         if (i.lt.k .or. k.eq.n) then
c
c     We cannot compute a reliable value for the gap of the last
c     Ritz value unless we know it is an approximation to the 
c     smallest singular value (k.eq.n). In this case we can take the 
c     distance to the next bigger one as the gap, which can really 
c     save us from getting stuck on matrices with a single isolated tiny 
c     singular value.
c
            if (i.eq.1) then
               gap = abs(theta(i)-theta(i+1))-max(bound(i),bound(i+1))
            else if (i.eq.n) then
               gap = abs(theta(i-1)-theta(i))-max(bound(i-1),bound(i))
            else
               gap = abs(theta(i)-theta(i+1))-max(bound(i),bound(i+1)) 
               gap = min(gap,abs(theta(i-1) - theta(i)) - 
     c              max(bound(i-1),bound(i)))
            endif
            if (gap.gt.bound(i)) then
               bound(i) = bound(i) * (bound(i)/gap)
            endif
         endif
      enddo
      end

      subroutine pdzero(n, x , incx)
      implicit none
      integer n, incx
      double precision x(*)
      
      integer i

      if ((n.gt.0).and.(incx.ne.0)) then         
         if (incx.eq.1) then
c$OMP PARALLEL DO 
            do i=1,n
               x(i) = 0d0
            enddo
         else
c$OMP PARALLEL DO firstprivate(incx) schedule(static)
            do i=1,n
               x(1+(i-1)*incx) = 0d0
            enddo
         endif
      endif
      return
      end

c     *****************************************************************
      subroutine pizero(n, x , incx)
      implicit none
      integer n, incx
      integer x(*)
      
      integer i

      if ((n.gt.0).and.(incx.ne.0)) then         
         if (incx.eq.1) then
c$OMP PARALLEL DO  schedule(static)
            do i=1,n
               x(i) = 0
            enddo
         else
c$OMP PARALLEL DO firstprivate(incx) schedule(static)
            do i=1,n
               x(1+(i-1)*incx) = 0
            enddo
         endif
      endif
      return
      end

c     *****************************************************************
      subroutine pdaxpy(n, alpha, x , incx, y, incy)
      implicit none
      integer n, incx, incy
      double precision alpha,x(*),y(*)
     
      call daxpy(n, alpha, x , incx, y, incy)
      end

c     *****************************************************************
      double precision function pddot(n, x , incx, y, incy)
      implicit none
      integer n, incx, incy
      double precision x(*),y(*),ddot
      external ddot
      
      pddot = ddot(n, x , incx, y, incy)
      end
      
c     *****************************************************************
      subroutine pdscal(n, alpha, x , incx)
      implicit none
      integer n, incx
      double precision alpha,x(*)

      call dscal(n, alpha, x , incx)
      end

c     *****************************************************************
      subroutine pdcopy(n, x , incx, y, incy)
      implicit none
      integer n, incx, incy
      double precision x(*),y(*)

      call dcopy(n, x , incx, y, incy)
      end