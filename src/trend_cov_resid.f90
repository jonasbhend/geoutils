SUBROUTINE trendcov(n, y, x, beta, betacov, missval, missthresh)

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
  !   betacov:    error covariance for slope of regression
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
  INTEGER ::              i, j, trndln, nspace
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: tmpval
  DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: xresid, tresid, yresid, temp  
  DOUBLE PRECISION, DIMENSION (n(1),n(2)) :: y, x
  DOUBLE PRECISION, DIMENSION (n(1)) :: beta
  DOUBLE PRECISION, DIMENSION (n(1),n(1)) :: betacov
  DOUBLE PRECISION :: missval, missthresh
  DOUBLE PRECISION :: sx, sy, ss, sxoss, a, b, st2, missnum
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: xt, yt, t

  ! Allocate temporary arrays (on the heap!)
  ALLOCATE(tmpval(n(1)))
  ALLOCATE(xt(n(2)))
  ALLOCATE(yt(n(2)))
  ALLOCATE(t(n(2)))
  ALLOCATE(temp(n(1),n(2)))
  ALLOCATE(yresid(n(1),n(2)))
  ALLOCATE(xresid(n(1),n(2)))
  ALLOCATE(tresid(n(1),n(2)))
 
  ! Initialization of values needed
  yresid(:,:) = 0
  xresid(:,:) = 0
  tresid(:,:) = 0
  temp(:,:) = 1
  trndln = n(2)
  nspace = n(1)
  beta(:) = missval
  betacov(:,:) = missval
  missnum = REAL(n(2)) * missthresh

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
        yresid(i,:) = yt - a - b*xt
        xresid(i,:) = t
     END IF
  END DO
  
  ! ! second pass to mask all missing values out
  ! t(:) = 1
  ! DO i = 1, nspace
  !    IF (tmpval(i) .GE. missnum .AND. tmpval(i) .LT. n(2)) THEN
  !       DO j = 1, nspace
  !          WHERE (y(i,:) .EQ. missval) 
  !             yresid(j,:) = 0
  !             xresid(j,:) = 0
  !          END WHERE
  !       END DO
  !       WHERE (y(i,:) .EQ. missval) t(:) = 0
  !   END IF
  ! END DO

  ! Last pass to compute correlation
  t(:) = 1
  DO i = 1, nspace
     IF (tmpval(i) .GE. missnum) THEN
        DO j = 1, i
           IF (tmpval(j) .GE. missnum) THEN
              betacov(i,j) = SUM(yresid(i,:)*yresid(j,:), MASK = y(i,:) .NE. missval .AND. y(j,:) .NE. missval ) / &
                   (SUM(xresid(i,:)*xresid(j,:), MASK = y(i,:) .NE. missval .AND. y(j,:) .NE. missval)) / &
                   (SUM(t*t, MASK = y(i,:) .NE. missval .AND. y(j,:) .NE. missval) - 1)
              betacov(j,i) = betacov(i,j)
           END IF
        END DO
     END IF
  END DO
  y = yresid
  x = xresid

  ! Deallocate temporary arrays
  DEALLOCATE(yresid)
  DEALLOCATE(xresid)
  DEALLOCATE(tresid)
  DEALLOCATE(t)
  DEALLOCATE(yt)
  DEALLOCATE(xt)

  DEALLOCATE(tmpval)

END SUBROUTINE trendcov
