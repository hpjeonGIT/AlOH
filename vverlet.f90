!############################################################################
SUBROUTINE VINIT(Nset, commn, param, angle, hbond, torsn, sys, q, cell, dt)
USE DATAFMT
IMPLICIT NONE
!
TYPE(AMNT):: Nset
TYPE(GNRL):: commn
TYPE(ATMC):: param(Natom)
TYPE(ANGL):: angle(Natom,Natom,Natom)
TYPE(HBND):: hbond(Natom)
TYPE(TRSN):: torsn(Natom,Natom,Natom,Natom)
TYPE(STAT):: sys
TYPE(PTCL):: q(Nset%Npt_all)
TYPE(SORT):: cell(Nset%Ncell_all)
REAL  (DP):: dt
!
TYPE(PTCL):: old(Nset%Npt_all)
TYPE(CORR):: pair
INTEGER(GI):: i, j, id, Npt(Natom)
REAL(DP):: xm, lambda(Natom), T0, xv(Nset%Npt_all,3), mv2(Natom), sumv(3), rand
!
! Velocity initialization by random number - rigid motion is removed
sumv(:) = 0.0D0
DO i = 1, Nset%Npt_all
   DO j=1, 3
      CALL RANDOM_NUMBER(rand)
      xv(i,j) = rand - 0.5D0
      sumv(j) = sumv(j) + xv(i,j)
   END DO
END DO
!
sumv(:) = sumv(:)/DBLE(Nset%Npt_all)
mv2(:) = 0.0D0
Npt(:) = 0
DO i=1, Nset%Npt_all
   xv(i,:)= xv(i,:) - sumv(:)
   id = q(i)%id
   xm = param(id)%am   
   mv2(id) = mv2(id) + xm*(xv(i,1)**2 + xv(i,2)**2 + xv(i,3)**2)
   Npt(id) = Npt(id) + 1
END DO
lambda(:) = 0.0D0
DO i=1, Natom
   IF (Npt(i) > 0) THEN
      T0 = mv2(i)/DBLE(Npt(i)*3)
      lambda(i) = DSQRT(sys%T_given/T0)
   END IF
END DO
!
DO i = 1, Nset%Npt_all
   id = q(i)%id
   q(i)%xv(:) = xv(i,:)*lambda(id)
   q(i)%ff(:) = 0.0D0
   old(i)%xr(:) = q(i)%xr(:)
   old(i)%xx(:) = q(i)%xx(:)
   old(i)%xv(:) = q(i)%xv(:)
END DO
!
! Initial force calculation - compensate the motion by the initial velocity
CALL VV_NVE1(Nset, param, q, -dt)
pair%dr = sys%Rcut/DBLE(Nrdf)
pair%rdf(:,:,:) = 0
CALL FORCE(Nset, commn, param, angle, hbond, torsn, sys, pair, cell, q)
sys%mv2 = 0.0D0
DO i=1, Nset%Npt_all
   id = q(i)%id
   xm = param(id)%am
   q(i)%xx(:) = old(i)%xx(:)
   q(i)%xv(:) = old(i)%xv(:)
   q(i)%xr(:) = old(i)%xr(:)
   sys%mv2 = sys%mv2 + xm*(q(i)%xv(1)**2 + q(i)%xv(2)**2 + q(i)%xv(3)**2 )
END DO
!
RETURN
END SUBROUTINE Vinit
!
!############################################################################
SUBROUTINE VV_NVE1(Nset, param, q, dt)
USE DATAFMT
!
TYPE(AMNT):: Nset
TYPE(ATMC):: param(Natom)
TYPE(PTCL):: q(Nset%Npt_all)
REAL(DP)  :: dt
!
INTEGER(GI):: i, id
REAL(DP):: xm, x(3), box(3)
!
box(:) = Nset%box(:)
DO i=1, Nset%Npt_all
   id = q(i)%id
   xm = param(id)%am
   q(i)%xv(:) = q(i)%xv(:) + 0.5D0 * dt * q(i)%ff(:)/xm
   x(:) = q(i)%xx(:) + dt*q(i)%xv(:)
   q(i)%xr(:) = q(i)%xr(:) + dt*q(i)%xv(:)
   q(i)%xx(:) = x(:) - box(:)*DNINT(x(:)/box(:))
END DO
!
RETURN
END SUBROUTINE VV_NVE1
!
!############################################################################
SUBROUTINE VV_NVE1_REFLECT(Nset, param, sys, q, dt)
USE DATAFMT
IMPLICIT NONE
!
TYPE(AMNT):: Nset
TYPE(ATMC):: param(Natom)
TYPE(STAT):: sys
TYPE(PTCL):: q(Nset%Npt_all)
REAL(DP)  :: dt
!
INTEGER(GI):: i, id
REAL(DP):: xm, x(3), box(3), box_z, tmp
!
box(:) = Nset%box(:)
box_z  = 0.5D0*box(3)
DO i=1, Nset%Npt_all
   id = q(i)%id
   xm = param(id)%am
   q(i)%xv(:) = q(i)%xv(:) + 0.5D0 * dt * q(i)%ff(:)/xm
   x(:) = q(i)%xx(:) + dt*q(i)%xv(:)
   q(i)%xx(1:2) = x(1:2) - box(1:2)*DNINT(x(1:2)/box(1:2))
   q(i)%xr(1:2) = q(i)%xr(1:2) + dt*q(i)%xv(1:2)
   !
   ! Reflecting boundary condition on the top
   IF (x(3) > box_z) THEN
      x(3) = box(3) - x(3)
      tmp = DSQRT(sys%T_given*3.D0/xm)
      q(i)%xv(3) = - tmp
      q(i)%xv(1:2) = 0.0D0
   END IF
   q(i)%xx(3) = x(3)
   q(i)%xr(3) = q(i)%xr(3) + dt*q(i)%xv(3)
END DO
!
RETURN
END SUBROUTINE VV_NVE1_REFLECT
!
!############################################################################
SUBROUTINE VV_NVE2(Nset, param, sys, q, dt)
USE DATAFMT
IMPLICIT NONE
!
TYPE(AMNT):: Nset
TYPE(ATMC):: param(Natom)
TYPE(STAT):: sys
TYPE(PTCL):: q(Nset%Npt_all)
REAL(DP)  :: dt
!
INTEGER(GI):: i, k, id
REAL(DP)   :: xm, mv2_tmp
!
sys%mv2 = 0.0D0
DO i=1, Nset%Npt_all
   id = q(i)%id
   xm = param(id)%am
   q(i)%xv(:) = q(i)%xv(:) + 0.5D0 * dt * q(i)%ff(:)/xm
   mv2_tmp = 0.0D0   
   DO k=1,3
      mv2_tmp = mv2_tmp + q(i)%xv(k)**2
   END DO
   sys%mv2 = sys%mv2 + xm*mv2_tmp
END DO
!
RETURN
END SUBROUTINE VV_NVE2
!
!############################################################################
SUBROUTINE VV_NVT2_BEREND(Nset, param, sys, q, dt)
USE DATAFMT
IMPLICIT NONE
!
TYPE(AMNT):: Nset
TYPE(ATMC):: param(Natom)
TYPE(STAT):: sys
TYPE(PTCL):: q(Nset%Npt_all)
REAL(DP)  :: dt
!
INTEGER(GI):: i, k, id, Npt(Natom)
REAL(DP):: xm, mv2_tmp, T0, lambda(Natom), mv2_sum(Natom)
REAL(DP), PARAMETER:: C_BEREND = 0.1
!
!
mv2_sum(:) = 0.0D0; lambda(:) = 0.0D0
Npt(:) = 0
DO i=1, Nset%Npt_all
   id = q(i)%id
   xm = param(id)%am
   q(i)%xv(:) = q(i)%xv(:) + 0.5D0 * dt * q(i)%ff(:)/xm
   mv2_tmp = 0.0D0
   DO k=1,3
      mv2_tmp = mv2_tmp + xm*q(i)%xv(k)**2
   END DO
   mv2_sum(id) = mv2_sum(id) + mv2_tmp
   Npt(id) = Npt(id) + 1
END DO
!
DO i=1, Natom
   IF (Npt(i)> 0) THEN
      T0 = mv2_sum(i)/DBLE(3*Npt(i))
      lambda(i) = DSQRT(1.D0 + C_BEREND*(sys%T_given/T0 - 1.0D0))
   END IF
END DO
!
sys%mv2 = 0.0D0
DO i=1, Nset%Npt_all
   id = q(i)%id
   xm = param(id)%am
   q(i)%xv(:) = q(i)%xv(:)*lambda(id)
   mv2_tmp = 0.0D0
   DO k=1,3
      mv2_tmp = mv2_tmp + xm*q(i)%xv(k)**2
   END DO
   sys%mv2 = sys%mv2 + mv2_tmp
END DO
!
RETURN
END SUBROUTINE VV_NVT2_BEREND
!
!############################################################################
SUBROUTINE VV_NVT1_STOCH(Nset, param, q, dt)
USE DATAFMT
IMPLICIT NONE
!
TYPE(AMNT):: Nset
TYPE(ATMC):: param(Natom)
TYPE(PTCL):: q(Nset%Npt_all)
REAL(DP)  :: dt
!
INTEGER(GI):: i, id
REAL(DP):: xm, x(3), box(3), alpha
!
alpha = Nset%alpha
box(:) = Nset%box(:)
DO i=1, Nset%Npt_all
   id = q(i)%id
   xm = param(id)%am
   q(i)%xv(:) = q(i)%xv(:)*(1.D0 - 0.5D0*alpha*dt/xm) + &
        & 0.5D0 * dt * q(i)%ff(:)/xm
   x(:) = q(i)%xx(:) + dt*q(i)%xv(:)
   q(i)%xr(:) = q(i)%xr(:) + dt*q(i)%xv(:)
   q(i)%xx(:) = x(:) - box(:)*DNINT(x(:)/box(:))
END DO
!
RETURN
END SUBROUTINE VV_NVT1_STOCH
!
!############################################################################
SUBROUTINE VV_NVT2_STOCH(Nset, param, sys, q, dt)
USE DATAFMT
IMPLICIT NONE
!
TYPE(AMNT):: Nset
TYPE(ATMC):: param(Natom)
TYPE(STAT):: sys
TYPE(PTCL):: q(Nset%Npt_all)
REAL(DP)  :: dt
!
INTEGER(GI):: i, k, id
REAL(DP):: xm, mv2_tmp, fluct, sigma, eta, alpha, beta
!
!
alpha = Nset%alpha
sigma = 1.0D0
beta = DSQRT(2.D0*alpha*sys%T_given/dt)
sys%mv2 = 0.0D0
DO i=1, Nset%Npt_all
   id = q(i)%id
   xm = param(id)%am
   mv2_tmp = 0.0D0
   DO k=1,3
      eta = fluct(sigma)
      q(i)%ff(k) = q(i)%ff(k) + eta*beta
      q(i)%xv(k) = (q(i)%xv(k) + 0.5D0*dt*q(i)%ff(k)/xm) / &
           &       (1.D0 + 0.5D0*alpha*dt/xm)
      mv2_tmp = mv2_tmp + xm*q(i)%xv(k)**2
   END DO
   sys%mv2 = sys%mv2 + mv2_tmp
END DO
!
!
RETURN
END SUBROUTINE VV_NVT2_STOCH
!
!############################################################################
SUBROUTINE REMOVE_RIGID_MOTION(Nset, q)
USE DATAFMT
IMPLICIT NONE
!
TYPE(AMNT):: Nset
TYPE(PTCL):: q(Nset%Npt_all)
!
INTEGER(GI):: i
REAL   (DP):: sumv(3)
!
sumv(:) = 0.0D0
DO i=1, Nset%Npt_all
   sumv(:) = sumv(:) + q(i)%xv(:)
END DO
sumv(:) = sumv(:)/DBLE(Nset%Npt_all)
DO i=1, Nset%Npt_all
   q(i)%xv(:) = q(i)%xv(:) - sumv(:)
END DO
!
RETURN
END SUBROUTINE REMOVE_RIGID_MOTION
