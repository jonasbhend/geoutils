SUBROUTINE trendresid(n, y, x, beta, dw, se0, seout, missval, missthresh)

  !   DESCRIPTION:
  !
  !   This subroutine computes trends of y on x and the respective covariance of
  !   the estimated slope using the residuals of the 
  !
  !   VARIABLES:
  !
  !   n:          vector with dimensions of y
  !   y:          predictand array, dim = ('space', time)
  !   x:          predictor array, dim = ('space', time)
  !   beta:       estimate of slope of regression
  !   missval:    missing value specification
  !   missthresh: threshold at which trends are computed (not more than
  !               1-missthresh values missing)
  !
  !   CREATED:    Oct. 23, 2006
  !
  !   MODIFIED:   July 11, 2011
  !
  !   OWNER:      Jonas Bhend, jonas.bhend_at_csiro.au
  !

  IMPLICIT NONE

  INTEGER (KIND = 4) ::   n(2)
  INTEGER ::              i, j, k, trndln, nspace
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: tmpval, dwcrit, acf
  DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: temp  
  DOUBLE PRECISION, DIMENSION (n(1),n(2)) :: y, x
  DOUBLE PRECISION, DIMENSION (n(1)) :: beta, dw, se0, seout
  DOUBLE PRECISION :: missval, missthresh
  DOUBLE PRECISION :: sx, sy, ss, sxoss, a, b, st2, missnum, sigb, sigdat
  DOUBLE PRECISION :: rho1, rho2, phi1, phi2, acsum, neff1, neff2, chi2
  DOUBLE PRECISION :: se1, se2
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: xt, yt, t, resid
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: resplus, resminus

  ! Allocate temporary arrays (on the heap!)
  ALLOCATE(tmpval(n(1)))
  ALLOCATE(dwcrit(7))
  ALLOCATE(acf(n(2)))

  ALLOCATE(xt(n(2)))
  ALLOCATE(yt(n(2)))
  ALLOCATE(t(n(2)))
  ALLOCATE(resid(n(2)))

  ALLOCATE(temp(n(1),n(2)))

  ALLOCATE(resplus(n(2) - 2))
  ALLOCATE(resminus(n(2) - 2))
 
  ! Initialization of values needed
  temp(:,:) = 1
  trndln = n(2)
  nspace = n(1)
  beta(:) = missval
  dw(:) = missval
  missnum = REAL(n(2)) * missthresh
  !   dwcrit = 0
  !   dwcrit(1) = 1.2
  !   dwcrit(2) = 1.29
  !   dwcrit(3) = 1.35
  !   dwcrit(4) = 1.44
  !   dwcrit(5) = 1.5
  !   dwcrit(6) = 1.55
  !   dwcrit(7) = 1.65
  

  ! calculating the slope of the regression
  DO i = 1,nspace
     WHERE (y(i,:) .EQ. missval) temp(i,:) = 0
     tmpval(i) = REAL(SUM(temp(i,:)))
     IF (tmpval(i) .GE. missnum) THEN
        yt = y(i,:) 
        xt = x(i,:)
        WHERE (y(i,:) .EQ. missval) &
             yt(:) = 0 
        WHERE (y(i,:) .EQ. missval) &
           xt(:) = 0 
        sx = SUM(xt)
        sy = SUM(yt)
        ss = INT(tmpval(i))
        sxoss = sx/ss
        t(:) = 0
        WHERE (y(i,:) .NE. missval) t(:) = xt(:) - sxoss
        st2 = SUM(t*t)
        b = SUM(t*yt)/st2
        a = (sy - sx*b)/ss
        beta(i) = b
        resid = yt - a - b*xt
        WHERE (y(i,:) .EQ. missval) resid(:) = 0
        chi2 = SUM(resid * resid)
        sigdat = sqrt(chi2 / (tmpval(i) - 2))              
        sigb = sqrt(1.0/st2) * sigdat
        se0(i) = sigb
       
         ! compute autocorrelation
        chi2 = SUM(resid*resid)
        rho1 = SUM(resid(2:trndln)*resid(1:(trndln-1)))/chi2
        resplus = resid(3:trndln) - SUM(resid(3:trndln))/(trndln - 2)
        resminus = resid(1:(trndln-2)) - SUM(resid(1:(trndln-2)))/(trndln - 2)
        rho2 = SUM(resplus*resminus)/sqrt(SUM(resplus*resplus)*SUM(resminus*resminus))
        
        ! AR(1) process
        acsum = 0
        DO k = 1, INT(tmpval(i)-1)
           acf(k) = rho1**k
           acsum  = acsum + (1 - k/tmpval(i))*acf(k)
        END DO
        neff1 = MIN(tmpval(i), tmpval(i)/(1 + 2*acsum))
        se1 = se0(i)*sqrt(tmpval(i))/sqrt(neff1)
        
       ! AR(2) process
        phi1 = rho1 * (1 - rho2) / (1 - rho1**2)
        phi2 = (rho2 - rho1**2) / (1 - rho1**2)
        acf(1) = rho1
        acf(2) = rho2
        acsum = (1 - 1/tmpval(i))*acf(1) + (1 - 2/tmpval(i))*acf(2)
        DO k = 3, INT(tmpval(i) - 1)
           acf(k) = phi1*acf(k-1) + phi2*acf(k-2)
           acsum = acsum + (1 - k/tmpval(i))*acf(k)
        END DO
        neff2 = MIN(tmpval(i), tmpval(i)/(1 + 2*acsum))
        se2 = se0(i) *sqrt(tmpval(i)) / sqrt(neff2)
        
       ! Decide which standard error to return
        seout(i) = se0(i)
        IF (rho1.GT.1.96/SQRT(tmpval(i))) THEN
           IF (ABS(phi2).GT.1.96/SQRT(tmpval(i))) THEN
              seout(i) = se2
           ELSE
              seout(i) = se1
           END IF
        END IF

        ! compute durbin-watson statistic
        dw(i) = SUM((resid(1:(trndln-1)) - resid(2:trndln))**2)/SUM(resid*resid)
        
        ! substitute missing values back in and write back to y
        WHERE (y(i,:) .EQ. missval) resid(:) = missval
        y(i,:) = resid
     ELSE 
        y(i,:) = missval   
     END IF
  END DO

  ! Deallocate temporary arrays
  DEALLOCATE(t)
  DEALLOCATE(yt)
  DEALLOCATE(resid)
  DEALLOCATE(xt)
  DEALLOCATE(acf)
  DEALLOCATE(dwcrit)
  DEALLOCATE(resminus)
  DEALLOCATE(resplus)

  DEALLOCATE(tmpval)

END SUBROUTINE trendresid
