SUBROUTINE artrend(n, y, x, trndln, yout, yresid, seout, se0, missval, missthresh, robust, gaps)

  !   DESCRIPTION:
  !
  !   This program runs a moving trend analysis based of cov/var
  !   for a single timeseries, Missing values are excluded
  !   it is invoked by movartrend in ~/R/geoUtils.R
  !   
  !   MODIFICATIONS:
  !   11/7/11     Test partial autocorrelation to decide what AR model to use
  !               PACF_crit ~= 2/sqrt(N)
  !
  !   VARIABLES:
  !
  !   n:          vector with dimensions of y
  !   y:          predictand array, dim = ('space', time)
  !   x:          predictor array, dim = ('space', time)
  !   trndln:     length of moving trend window (normally uneven,
  !               trends assigned to LAST timestep)
  !   yout:       regression slope estimates, dim(y)
  !   yresid:     residuals from regression (last if movtrend)
  !   seout:      standard error for slope, dim(y)
  !   se0:        standard error, no autocorrelation
  !   se1:        standard error corrected for effective sample size AR(1), dim(y)
  !   se2:        standard error corrected for effective sample size AR(2), dim(y)
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
  INTEGER (KIND = 4) ::   trndln
  INTEGER ::              i, j, k, l, c, c2, c3, ntmp
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: mean, tmpval
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: ymean1, ymean2
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: tmpval1, tmpval2
  DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: tmpgap
  DOUBLE PRECISION, DIMENSION (n(1),n(2)) :: y, x, yout, yresid, seout, se0
  DOUBLE PRECISION :: se1, se2, neff1, neff2
  DOUBLE PRECISION :: missval
  DOUBLE PRECISION :: sx, sy, ss, sxoss, a, b, sigb, sigdat, chi2, st2
  DOUBLE PRECISION :: missthresh, missthresh1, missthresh2
  DOUBLE PRECISION :: phi1, phi2, rho1, rho2, acsum
  DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: ytemp, xtemp, temp
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: xt, yt, resid, acf
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: resplus, resminus
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: t
  INTEGER, DIMENSION (:), ALLOCATABLE :: index
  INTEGER (KIND = 4):: robust, gaps
  REAL, PARAMETER :: TINY = 1e-20

  ! Allocate temporary arrays (on the heap!)
  ALLOCATE(mean(n(1)))
  ALLOCATE(tmpval(n(1)))

  ALLOCATE(ymean1(n(1)))
  ALLOCATE(ymean2(n(1)))
  
  ALLOCATE(tmpval1(n(1)))
  ALLOCATE(tmpval2(n(1)))

  ALLOCATE(tmpgap(n(1),(trndln-1)))

  ALLOCATE(ytemp(n(1),trndln))
  ALLOCATE(xtemp(n(1),trndln))
  ALLOCATE(temp(n(1),trndln))

  ALLOCATE(xt(trndln))
  ALLOCATE(yt(trndln))
  ALLOCATE(resid(trndln))
  ALLOCATE(acf(trndln))

  ALLOCATE(resplus(trndln - 2))
  ALLOCATE(resminus(trndln - 2))

  ALLOCATE(t(trndln))
  
  ALLOCATE(index(n(1)))

  ! Initialization of values needed
  temp = 1
  c = trndln - 1
  yout = missval
  seout = missval
  neff1 = missval
  neff2 = missval
  tmpgap = 0   


  SELECT CASE (robust)

  CASE (1)

     c2 = INT(c/2)
     c3 = INT(trndln/2)
     missthresh1 = missthresh * REAL(c3)
     missthresh2 = missthresh * REAL(c2 +1) 

     DO j = trndln,n(2)
        ytemp(:,1:c3) = y(:,(j-c):(j-trndln +c3))
        ymean1 = SUM(ytemp(:,1:c3), DIM=2, MASK=ytemp(:,1:c3) .NE. missval)
        tmpval1 = SUM( temp(:,1:c3), DIM=2, MASK=ytemp(:,1:c3) .NE. missval )
        ymean1 = ymean1 / (tmpval1 + TINY)

        ytemp(:,1:(c2+1)) = y(:,(j-c2):j)
        ymean2 = SUM(ytemp(:,1:(c2+1)), DIM=2, &
             MASK=ytemp(:,1:(c2+1)) .NE. missval)
        tmpval2 = SUM( temp(:,1:(c2+1)), DIM=2, &
             MASK=ytemp(:,1:(c2+1)) .NE. missval )
        ymean2 = ymean2 / (tmpval2 + TINY)

        yout(:,j) = (ymean2 - ymean1) / (REAL(trndln) /2)
        WHERE (tmpval1(:) .LT. missthresh1 .OR. tmpval2(:) .LT. missthresh2) &
             yout(:,j) = missval

     ENDDO

  CASE (0)

     missthresh = REAL(trndln) * missthresh

     DO j = trndln,n(2)
        ytemp = y(:,(j-c):j)
        xtemp = x(:,(j-c):j)
        tmpval = SUM( temp, DIM=2, MASK=ytemp(:,:) .NE. missval )

        ! select which i to loop on, write to index(1:ntmp)
        IF (gaps.EQ.1) THEN
           tmpgap(:,:) = 0

           WHERE (ytemp(:,1:c) .EQ. missval) tmpgap(:,:) = 1
           WHERE (ytemp(:,2:trndln) .EQ. missval) tmpgap(:,:) = tmpgap(:,:) + 1

           WHERE (ANY(tmpgap(:,:) .GT. 1, DIM = 2)) tmpval(:) = 0

        END IF


        ! calculating the slope of the regression
        DO i = 1, n(1)
           IF (tmpval(i) .GT. missthresh) THEN
              yt = ytemp(i,:) 
              xt = xtemp(i,:)
              WHERE (ytemp(i,:) .EQ. missval .OR. xtemp(i,:) .EQ. missval) &
                   yt(:) = 0 
              WHERE (ytemp(i,:) .EQ. missval .OR. xtemp(i,:) .EQ. missval) &
                   xt(:) = 0 
              sx = SUM(xt)
              sy = SUM(yt)
              ss = INT(tmpval(i))
              sxoss = sx/ss
              t(:) = 0
              WHERE (ytemp(i,:) .NE. missval .AND. xtemp(i,:) .NE. missval) t(:) = xt(:) - sxoss
              st2 = SUM(t*t)
              b = SUM(t*yt)/st2
              a = (sy - sx*b)/ss
              resid(:) = 0
              WHERE (ytemp(i,:) .NE. missval .AND. xtemp(i,:) .NE. missval) resid(:) = yt(:) - a - b*xt(:)
              chi2 = SUM(resid*resid)
              sigdat = sqrt(chi2 / (tmpval(i) - 2))              
              sigb = sqrt(1.0/st2) * sigdat
              
              yout(i,j) = b
              yresid(i,(j - trndln + 1):j) = resid
              WHERE (ytemp(i,:) .EQ. missval .OR. xtemp(i,:) .EQ. missval) yresid(i,(j - trndln + 1):j) = missval
              se0(i,j) = sigb

              ! compute autocorrelation
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
              ! neff1 = MIN(tmpval(i), tmpval(i)/(1 + 2*acsum))
              ! se1 = se0(i,j)*sqrt(tmpval(i))/sqrt(neff1)
              se1 = se0(i,j) * sqrt(1 + 2*acsum)
              
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
              ! neff2 = MIN(tmpval(i), tmpval(i)/(1 + 2*acsum))
              ! se2 = se0(i,j) *sqrt(tmpval(i)) / sqrt(neff2)
              se2 = se0(i,j) * sqrt(1 + 2*acsum)

              seout(i,j) = se1
              
              ! test for autocorrelation
              !IF (ABS(rho1).GT.2/SQRT(tmpval(i))) THEN
              !   IF (ABS(phi2).GT.2/SQRT(tmpval(i))) THEN
              !      seout(i,j) = se2
              !   ELSE
              !      seout(i,j) = se1
              !   END IF
              !END IF
              
           END IF
           
        ENDDO

     ENDDO

  END SELECT

  ! Deallocate temporary arrays
  DEALLOCATE(index)
  
  DEALLOCATE(t)

  DEALLOCATE(resminus)
  DEALLOCATE(resplus)

  DEALLOCATE(acf)
  DEALLOCATE(resid)
  DEALLOCATE(yt)
  DEALLOCATE(xt)

  DEALLOCATE(temp)
  DEALLOCATE(xtemp)
  DEALLOCATE(ytemp)

  DEALLOCATE(tmpgap)

  DEALLOCATE(tmpval2)
  DEALLOCATE(tmpval1)
  
  DEALLOCATE(ymean2)
  DEALLOCATE(ymean1)

  DEALLOCATE(tmpval)
  DEALLOCATE(mean)

END SUBROUTINE artrend
