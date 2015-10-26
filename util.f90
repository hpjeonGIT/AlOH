!
!#############################################################################
SUBROUTINE PRINT_STAT(Nset, sys, Nloop, cell, q, dt)
USE DATAFMT
IMPLICIT NONE
TYPE(AMNT):: Nset
TYPE(STAT):: sys
TYPE(SORT):: cell(Nset%Ncell_all)
TYPE(PTCL):: q(Nset%Npt_all)
INTEGER(GI):: Nloop
REAL(DP):: dt
!
INTEGER(GI):: i, j, k, Npt_local, ni, nj, id_pair, m, n1, Npt_pair, id
REAL   (DP):: E_total, Temperature, P_short, dx(3), df(3), r2, &
     &        rcut2, rx(3), box(3)
REAL   (DP), PARAMETER:: kcalpmol = 23.06035D0
!
!
E_total = sys%Ereaxff + sys%mv2/2.D0
Temperature = sys%mv2/DBLE(Nset%Npt_all*3)/kB
rcut2 = sys%Rcut2
box(1:3) = sys%box(1:3)
P_short = 0.0D0
!
! Pressure calculations
!
!
WRITE(15,100) Nloop*TFM*dt, Temperature, E_total, sys%Ereaxff, sys%press
100 FORMAT (ES11.4, 1X, F6.1, 3(1X, ES11.4))
WRITE(17,200) Nloop*TFM*dt, sys%E_bond, sys%E_over, sys%E_under,  sys%E_lp, &
     & sys%E_val,  sys%E_pen,  sys%E_H, sys%E_tor, sys%E_conj,  sys%E_vdw, &
     & sys%E_coul, sys%E_charge*kcalpmol
200 FORMAT(13(ES10.3, 1X))
RETURN
END SUBROUTINE PRINT_STAT
!
!#############################################################################
SUBROUTINE PRINT_SNAPSHOT(Nset, Nloop, q)
USE DATAFMT
IMPLICIT NONE
TYPE(AMNT):: Nset
TYPE(PTCL):: q(Nset%Npt_all)
INTEGER(GI):: Nloop
!
INTEGER(GI):: openstatus, nframe, i, id
CHARACTER(LEN=256):: dummy, RESTFILE, SNAPFILE
!
nframe = Nloop/Nset%freq_snap
WRITE(RESTFILE,100) nframe
100 FORMAT("snapshot", I3.3, ".xyz")
WRITE(SNAPFILE,110) nframe
110 FORMAT("image", I3.3, ".xyz")
!
OPEN(UNIT=20, FILE=RESTFILE, ACTION="WRITE", &
     & FORM = "FORMATTED", IOSTAT = openstatus)
IF (openstatus >0) THEN
   PRINT '("=== cannot open ", A, "===")', RESTFILE ; STOP
END IF
OPEN(UNIT=21, FILE=SNAPFILE, ACTION="WRITE", &
     & FORM = "FORMATTED", IOSTAT = openstatus)
IF (openstatus >0) THEN
   PRINT '("=== cannot open ", A, "===")', SNAPFILE ; STOP
END IF

WRITE(20,*) Nset%Npt_all
WRITE(20,*) "Frame = ", Nframe, " energy =  0.0 "
WRITE(21,*) Nset%Npt_all
WRITE(21,*) "Frame = ", Nframe, " energy =  0.0 "

DO i=1, Nset%Npt_all
   id = q(i)%id
   SELECT CASE (id)
      CASE (ID_H)
         dummy = "H "
      CASE (ID_O)
         dummy = "O "
      CASE (ID_Al)
         dummy = "Al "
      CASE DEFAULT
         CLOSE(20)
         PRINT *, id
         STOP "=== Particle id error ==="
   END SELECT
   IF (Nset%tag_msd) THEN
      WRITE(20,200) dummy, q(i)%xr(:), q(i)%xv(:), q(i)%ff(:), q(i)%q
   ELSE
      WRITE(20,200) dummy, q(i)%xx(:), q(i)%xv(:), q(i)%ff(:), q(i)%q
   END IF
   WRITE(21,300) dummy, q(i)%xx(:), q(i)%q
END DO
200 FORMAT(A3, 3(ES13.6, 1X), 6(ES11.4, 1X), ES13.6)
300 FORMAT(A3, 4(ES11.4, 1X))
CLOSE(20)
CLOSE(21)
RETURN
END SUBROUTINE PRINT_SNAPSHOT
!
!****************************************************************************** 
FUNCTION fluct(x)
  USE DATAFMT
  IMPLICIT NONE
  REAL(DP):: fluct, x, r, v1, v2
  REAL(DP):: rand1, rand2
  !
  ! Initialization
  r=1.D0
  DO WHILE (r.ge.1.D0)
     CALL RANDOM_NUMBER(rand1)
     CALL RANDOM_NUMBER(rand2)
     v1 = 2.D0*rand1 - 1.D0
     v2 = 2.D0*rand2 - 1.D0
     r = v1*v1+v2*v2
  END DO
  fluct = v1*DSQRT(-2.D0*DLOG(r)/r)*x
  RETURN
END FUNCTION fluct
!
!###############################################################################
SUBROUTINE RANDOM_POSITION_O2(z_new, Lx, Ly, xx)
USE DATAFMT, ONLY: DP, GI, PI
REAL(DP):: z_new, Lx, Ly, xx(6)
!
REAL(DP):: cx(3), rand, phi, theta, L0
CALL RANDOM_NUMBER(rand)
cx(1) = (rand - 0.5D0)*Lx
CALL RANDOM_NUMBER(rand)
cx(2) = (rand - 0.5D0)*Ly
cx(3) = z_new
CALL RANDOM_NUMBER(rand)
phi = rand*PI*2.D0
L0 = 3.5D0
xx(1) = L0*DCOS(phi)
xx(2) = L0*DSIN(phi)
xx(3) = 0.0D0
xx(4:6) = - xx(1:3)
xx(1:3) = xx(1:3) + cx(1:3)
xx(1) = xx(1) - Lx*DNINT(xx(1)/Lx)
xx(2) = xx(2) - Ly*DNINT(xx(2)/Ly)
xx(4:6) = xx(4:6) + cx(1:3)
xx(4) = xx(4) - Lx*DNINT(xx(4)/Lx)
xx(5) = xx(5) - Ly*DNINT(xx(5)/Ly)
!
RETURN
END SUBROUTINE RANDOM_POSITION_O2
!
!###############################################################################
SUBROUTINE RANDOM_POSITION_OH(z_new, Lx, Ly, xx)
USE DATAFMT, ONLY: DP, GI, PI
REAL(DP):: z_new, Lx, Ly, xx(6)
!
REAL(DP):: cx(3), rand, phi, theta, L0
CALL RANDOM_NUMBER(rand)
cx(1) = (rand - 0.5D0)*Lx
CALL RANDOM_NUMBER(rand)
cx(2) = (rand - 0.5D0)*Ly
cx(3) = z_new
CALL RANDOM_NUMBER(rand)
phi = rand*PI*2.D0
L0 = 0.5D0
xx(1) = L0*DCOS(phi)
xx(2) = L0*DSIN(phi)
xx(3) = 0.0D0
xx(4:6) = - xx(1:3)
xx(1:3) = xx(1:3) + cx(1:3)
xx(1) = xx(1) - Lx*DNINT(xx(1)/Lx)
xx(2) = xx(2) - Ly*DNINT(xx(2)/Ly)
xx(4:6) = xx(4:6) + cx(1:3)
xx(4) = xx(4) - Lx*DNINT(xx(4)/Lx)
xx(5) = xx(5) - Ly*DNINT(xx(5)/Ly)
!
RETURN
END SUBROUTINE RANDOM_POSITION_OH
!
!###############################################################################
SUBROUTINE RANDOM_POSITION_O4(z_new, Lx, Ly, xx)
USE DATAFMT, ONLY: DP, GI, PI
REAL(DP):: z_new, Lx, Ly, xx(12)
!
REAL(DP):: cx(3), rand
CALL RANDOM_NUMBER(rand)
cx(1) = (rand - 0.5D0)*Lx
CALL RANDOM_NUMBER(rand)
cx(2) = (rand - 0.5D0)*Ly
cx(3) = z_new
L0 = 4.0D0
xx( 1) = cx(1) + L0; xx( 2) = cx(2) + L0; xx( 3) = cx(3)
xx( 4) = cx(1) + L0; xx( 5) = cx(2) - L0; xx( 6) = cx(3)
xx( 7) = cx(1) - L0; xx( 8) = cx(2) + L0; xx( 9) = cx(3)
xx(10) = cx(1) - L0; xx(11) = cx(2) - L0; xx(12) = cx(3)
xx( 1) = xx( 1) - Lx*DNINT(xx( 1)/Lx); xx( 2) = xx( 2) - Ly*DNINT(xx( 2)/Ly)
xx( 4) = xx( 4) - Lx*DNINT(xx( 4)/Lx); xx( 5) = xx( 5) - Ly*DNINT(xx( 5)/Ly)
xx( 7) = xx( 7) - Lx*DNINT(xx( 7)/Lx); xx( 8) = xx( 8) - Ly*DNINT(xx( 8)/Ly)
xx(10) = xx(10) - Lx*DNINT(xx(10)/Lx); xx(11) = xx(11) - Ly*DNINT(xx(11)/Ly)
!
RETURN
END SUBROUTINE RANDOM_POSITION_O4
!
!#############################################################################
SUBROUTINE OXIDE_CHECK(Nset, tag, q)
USE DATAFMT
IMPLICIT NONE
TYPE(AMNT):: Nset
TYPE(TCKT):: tag
TYPE(PTCL):: q(Nset%Npt_all)
!
INTEGER(GI):: i, ID_check
REAL   (DP):: o_height, q_max
!
! Highest metal position
Nset%al_zmax = 0.0D0
DO i=1, Nset%Npt_all
   IF (q(i)%id == ID_Al) THEN 
      Nset%al_zmax = MAX(q(i)%xx(3), Nset%al_zmax)
   END IF
END DO
ID_check = ID_O
!
! Only newly inserted O2
o_height = 0.0D0
q_max = -10.0D0
DO i=Nset%Npt_all, Nset%Npt_all-1, -1
   IF (q(i)%id == ID_check) THEN 
      o_height = MAX(q(i)%xx(3), o_height)
      q_max    = MAX(q(i)%q, q_max)
   END IF
END DO
!
IF ((Nset%al_zmax + 2.0D0) > o_height .OR. q_max < Nset%q_crit) THEN
   tag%new_o2 = .TRUE.
ELSE
   tag%new_o2 = .FALSE.
END IF
RETURN
END SUBROUTINE OXIDE_CHECK
!
!#############################################################################
SUBROUTINE PRINT_RDF(Nset, Nloop, pair)
USE DATAFMT
IMPLICIT NONE
TYPE(AMNT):: Nset
INTEGER(GI):: Nloop
TYPE(CORR):: pair
!
INTEGER(GI):: openstatus, nframe, i, j, n, Nkind, Npt_old, Npt_new
REAL   (DP):: result(Natom*(Natom+1)/2, Nrdf), r, dr, const, loop_sum, Xrdf
CHARACTER(LEN=256):: FILENAME  

!
nframe = Nloop/Nset%freq_rdf
WRITE(FILENAME,100) nframe
100 FORMAT("rdf", I3.3, ".dat")
!
OPEN(UNIT=20, FILE=FILENAME, ACTION="WRITE", &
     & FORM = "FORMATTED", IOSTAT = openstatus)
IF (openstatus >0) THEN
   PRINT '("=== cannot open ", A, "===")', FILENAME ; STOP
END IF
!
Xrdf = DBLE(Nset%freq_rdf)
n = 0
DO i=1, Natom
   DO j=i, Natom
      n = n + 1
      IF (i /= j) THEN
         result(n, 1:Nrdf) = DBLE(pair%rdf(i,j,1:Nrdf) + pair%rdf(j,i,1:Nrdf)) &
              &  /Xrdf
      ELSE
         result(n, 1:Nrdf) = DBLE(pair%rdf(i,j,1:Nrdf))/Xrdf
      END IF
   END DO
END DO
dr = pair%dr
Nkind = Natom*(Natom+1)/2
Npt_old = pair%Npt_old
Npt_new = pair%Npt_new
loop_sum = DBLE(Npt_old*Npt_new)*4.D0*PI/3.0D0/pair%V
WRITE(20,300) Nloop, Npt_old, Npt_new
300 FORMAT("# Header: current loop = ",I7," Npt_old = ", I5, " Npt_new = ", I5)
WRITE(20,400)
400 FORMAT("# Radius(A) H-H (2) / H-O (3) / H-Al (4) / O-O (5) / O-Al (6) " &
         & "/ Al-Al (7) ")
DO i=1, Nrdf
   r = DBLE(i)*dr
   const = ((r+dr)**3 - dr**3)*loop_sum
   WRITE(20,500, ADVANCE="NO") dr*(DBLE(i)-0.5D0)
   DO j=1, Nkind
      WRITE(20,500, ADVANCE="NO") result(j,i)/const
   END DO
   WRITE(20,*)
END DO
500 FORMAT(ES10.3, 1X)
CLOSE(20)
pair%rdf(:,:,:) = 0
pair%Npt_old = pair%Npt_new
!
RETURN
END SUBROUTINE PRINT_RDF
!
!###############################################################################
SUBROUTINE POTENTIAL_SETUP(commn, param, angle, hbond, torsn, sys)
USE DATAFMT
IMPLICIT NONE
TYPE(GNRL):: commn
TYPE(ATMC):: param(Natom)
TYPE(ANGL):: angle(Natom,Natom,Natom)
TYPE(HBND):: hbond(Natom)
TYPE(TRSN):: torsn(Natom,Natom,Natom,Natom)
TYPE (STAT):: sys
!
INTEGER   (GI):: openstatus, Nitem, i, j, k, m, n, n1, n2
CHARACTER(256):: dummy, name(Natom)
REAL      (DP):: Ra, Rb, rsig, rpi, rpipi, v1, v2, v3, v2BO, vconj, &
     &           ThetaO, kaa, kbb, pconj, pv2, kpenal, pv3
!REAL      (DP), PARAMETER:: kcalpmol = 23.06035D0
!
OPEN(UNIT=15,FILE="reaxff.aloh", STATUS="OLD", ACTION="READ", &
     & FORM = "FORMATTED", IOSTAT = openstatus)
IF (openstatus >0) STOP "### Cannot open reaxff.aloh ###"
READ(15,*) dummy
READ(15,*) dummy
!
! General parameter parsing
READ(15,*) Nitem
IF (Nitem /= 39) STOP "### Number of general parameter has been changed ###"
READ(15,*) commn%over1  ! #01
READ(15,*) commn%over2
READ(15,*) commn%val_ang1 
READ(15,*) commn%tri_bond1 
READ(15,*) commn%tri_bond2 
READ(15,*) commn%c2_corr 
READ(15,*) commn%under1 
READ(15,*) commn%tri_bond3 
READ(15,*) commn%under2 
READ(15,*) commn%under3
READ(15,*) commn%tri_bond4 !#11
READ(15,*) commn%l_taper  
READ(15,*) commn%u_taper 
READ(15,*) dummy
READ(15,*) commn%val_under 
READ(15,*) commn%val_anlp 
READ(15,*) commn%val_ang2 
READ(15,*) commn%val_ang3 
READ(15,*) dummy
READ(15,*) commn%dbl_bond1 
READ(15,*) commn%dbl_bond2 !#21
READ(15,*) commn%dbl_bond3 
READ(15,*) dummy
READ(15,*) commn%tor_BO 
READ(15,*) commn%tor_over1 
READ(15,*) commn%tor_over2 
READ(15,*) commn%conjug0 
READ(15,*) commn%conjug1 
READ(15,*) commn%vdw_sh 
commn%inv_vdwsh = 1.D0/commn%vdw_sh
READ(15,*) commn%cutoff 
READ(15,*) commn%val_ang_con1 !#31 
READ(15,*) commn%over3 
READ(15,*) commn%over4 
READ(15,*) commn%val_lp 
READ(15,*) dummy
READ(15,*) dummy
READ(15,*) dummy
READ(15,*) dummy
READ(15,*) commn%val_ang_con2 
commn%rcut_bo = commn%cutoff*0.01D0
!
! Taper correction coefficients
Ra = commn%l_taper  
Rb =  commn%u_taper 
sys%Rcut = Rb
commn%Taper(1) = -35.D0*Ra**3*Rb**4 + 21.D0*Ra**2*Rb**5 + 7.D0*Ra*Rb**6 + Rb**7
commn%Taper(2) = 140.D0*Ra**3*Rb**3
commn%Taper(3) = -210.D0*(Ra**3*Rb**2 + Ra**2*Rb**3)
commn%Taper(4) = 140.D0*(Ra**3*Rb + 3.D0*Ra**2*Rb**2 + Ra*Rb**3)
commn%Taper(5) = -35.D0*(Ra**3 + 9.D0*Ra**2*Rb + 9.D0*Ra*Rb**2 + Rb**3)
commn%Taper(6) = 84.D0*(Ra**2 + 3.D0*Ra*Rb + Rb**2)
commn%Taper(7) = -70.D0*(Ra + Rb)
commn%Taper(8) = 20.D0
commn%Taper(:) = commn%Taper(:)/(Rb - Ra)**7
!
! Atom data parsing
READ(15,*) Nitem
IF (Nitem /= Natom) STOP "### Number of atoms does not match ###"
READ(15,*) dummy;READ(15,*) dummy;READ(15,*) dummy;
DO i=1, Natom
   READ(15,*) name(i), param(i)%cov_r, param(i)%valency, param(i)%am, &
        & param(i)%Rvdw, param(i)%Evdw, param(i)%gammaEEM, param(i)%cov_r2, &
        & param(i)%no_el
   READ(15,*) param(i)%alfa, param(i)%gammavdw, param(i)%valency2, &
        & param(i)%Eunder, dummy, param(i)%chiEEM, param(i)%etaEEM
   READ(15,*) param(i)%cov_r3, param(i)%Elp, param(i)%Hinc, &
        & param(i)%BO13_1, param(i)%BO13_2, param(i)%BO13_3
   READ(15,*) param(i)%ovun, param(i)%val1, dummy, param(i)%val3, param(i)%vval4
   param(i)%nlp_opt = (param(i)%no_el - param(i)%valency)/2.D0
END DO
!
IF (name(1) /= "H" .OR. name(2) /= "O" .OR. name(3) /= "Al") &
     &  STOP "=== Error in element type ==="
DO i=1, Natom
   DO j=1, Natom
      commn%Rvdw(i,j)  = DSQRT(4.D0*param(i)%Rvdw*param(j)%Rvdw) !<-- check x2 or not
      commn%Dvdw(i,j)  = DSQRT(param(i)%Evdw*param(j)%Evdw)
      commn%alpha(i,j) = DSQRT(param(i)%alfa*param(j)%alfa)
      IF (param(i)%cov_r > ZZERO .AND. param(j)%cov_r > ZZERO) THEN
         commn%cov_r(i,j)= (param(i)%cov_r  + param(j)%cov_r )/2.0D0
      ELSE
         commn%cov_r(i,j) = -1.0
      END IF
      IF (param(i)%cov_r2 > ZZERO .AND. param(j)%cov_r2 > ZZERO) THEN
         commn%cov_r2(i,j)= (param(i)%cov_r2 + param(j)%cov_r2)/2.0D0
      ELSE
         commn%cov_r2(i,j) = -1.0
      END IF
      IF (param(i)%cov_r3 > ZZERO .AND. param(j)%cov_r3 > ZZERO) THEN
         commn%cov_r3(i,j)= (param(i)%cov_r3 + param(j)%cov_r3)/2.0D0
      ELSE
         commn%cov_r3(i,j) = -1.0
      END IF      
      commn%g_vdw(i,j) = (1.D0/DSQRT(param(i)%gammavdw* &
           & param(j)%gammavdw))**commn%vdw_sh
      commn%g_coul(i,j) = 1.D0/(DSQRT(param(i)%gammaEEM*param(j)%gammaEEM))**3
      commn%BO13_1(i,j) = DSQRT(param(i)%BO13_1*param(j)%BO13_1)
      commn%BO13_2(i,j) = DSQRT(param(i)%BO13_2*param(j)%BO13_2)
      commn%BO13_3(i,j) = DSQRT(param(i)%BO13_3*param(j)%BO13_3)
   END DO
END DO
!
! Bond data parsing
DO i=1, Natom
   param(i)%Edis1(:) = 0.0D0 ! Criterion for the bond-calculation
END DO
READ(15,*) Nitem
!IF (Nitem /= (Natom+1)*Natom/2) STOP "### bond data does not match ###"
READ(15,*) dummy
DO i=1, Nitem
   READ(15,*) m, n, param(m)%Edis1(n), param(m)%Edis2(n), param(m)%Edis3(n), &
        & param(m)%pbe1(n), param(m)%pbo5(n), param(m)%corr13(n), &
        & param(m)%pbo6(n), param(m)%kov(n)
   READ(15,*) param(m)%pbe2(n), param(m)%pbo3(n), param(m)%pbo4(n), dummy, &
        & param(m)%pbo1(n), param(m)%pbo2(n), param(m)%ovcorr(n)
   IF (m /= n) THEN
      param(n)%Edis1(m)  = param(m)%Edis1(n)
      param(n)%Edis2(m)  = param(m)%Edis2(n)
      param(n)%Edis3(m)  = param(m)%Edis3(n)
      param(n)%pbe1(m)   = param(m)%pbe1(n)
      param(n)%pbo5(m)   = param(m)%pbo5(n)
      param(n)%corr13(m) = param(m)%corr13(n)
      param(n)%pbo6(m)   = param(m)%pbo6(n)
      param(n)%kov(m)    = param(m)%kov(n)
      param(n)%pbe2(m)   = param(m)%pbe2(n)
      param(n)%pbo3(m)   = param(m)%pbo3(n)
      param(n)%pbo4(m)   = param(m)%pbo4(n)
      param(n)%pbo1(m)   = param(m)%pbo1(n)
      param(n)%pbo2(m)   = param(m)%pbo2(n)
      param(n)%ovcorr(m) = param(m)%ovcorr(n)
   END IF
END DO
!
! Offdiagonal data (exception for linear combination) parsing
READ(15,*) Nitem
DO i=1, Nitem
   READ(15,*) m, n, commn%Dvdw(m,n), commn%Rvdw(m,n), commn%alpha(m,n), &
        & rsig, rpi, rpipi
   IF (rsig > ZZERO) commn%cov_r(m,n) = rsig
   IF (rpi  > ZZERO) commn%cov_r2(m,n) = rpi
   IF (rpipi> ZZERO) commn%cov_r3(m,n) = rpipi
   commn%Rvdw(m,n)  = commn%Rvdw(m,n)*2.D0
   !
   commn%Dvdw(n,m)  = commn%Dvdw(m,n)
   commn%Rvdw(n,m)  = commn%Rvdw(m,n)
   commn%alpha(n,m) = commn%alpha(m,n)
   commn%cov_r(n,m) = commn%cov_r(m,n)
   commn%cov_r2(n,m)= commn%cov_r2(m,n)
   commn%cov_r3(n,m)= commn%cov_r3(m,n)
END DO
!
! Angle potential data parsing
DO i=1, Nang
   angle(:,:,:)%ThetaO(i) = 0.0D0
   angle(:,:,:)%ka    (i) = 0.0D0 ! Criterion for the angle-calculation
   angle(:,:,:)%kb    (i) = 0.0D0
   angle(:,:,:)%pconj (i) = 0.0D0
   angle(:,:,:)%pv2   (i) = 0.0D0
   angle(:,:,:)%kpenal(i) = 0.0D0
   angle(:,:,:)%pv3   (i) = 0.0D0
END DO
angle(:,:,:)%Nset = 0
!
READ(15,*) Nitem
DO i=1, Nitem
   READ(15,*) m, j, n, ThetaO, kaa, kbb, pconj, pv2, kpenal, pv3
   k = angle(m,j,n)%Nset 
   k = k+1
   IF (k > Nang) STOP "=== Number valency angle components crashed ==="
   ThetaO = ThetaO*PI/180.D0
   angle(m,j,n)%Nset = k
   angle(m,j,n)%Thetao(k) = ThetaO
   angle(m,j,n)%ka(k)     = kaa
   angle(m,j,n)%kb(k)     = kbb
   angle(m,j,n)%pconj(k)  = pconj
   angle(m,j,n)%pv2(k)    = pv2
   angle(m,j,n)%kpenal(k) = kpenal
   angle(m,j,n)%pv3(k)    = pv3
   IF (m /= n) THEN
      angle(n,j,m)%Nset = k
      angle(n,j,m)%ThetaO(k) = angle(m,j,n)%ThetaO(k)
      angle(n,j,m)%ka    (k) = angle(m,j,n)%ka    (k)
      angle(n,j,m)%kb    (k) = angle(m,j,n)%kb    (k)
      angle(n,j,m)%pconj (k) = angle(m,j,n)%pconj (k)
      angle(n,j,m)%pv2   (k) = angle(m,j,n)%pv2   (k)
      angle(n,j,m)%kpenal(k) = angle(m,j,n)%kpenal(k)
      angle(n,j,m)%pv3   (k) = angle(m,j,n)%pv3   (k)
   END IF
END DO
!
! Torsion potential data parsing
torsn(:,:,:,:)%V1 = 0.0D0
torsn(:,:,:,:)%V2 = 0.0D0
torsn(:,:,:,:)%V3 = 0.0D0  ! Criterion for the torsion-calculation
torsn(:,:,:,:)%V2_BO = 0.0D0
torsn(:,:,:,:)%vconj = 0.0D0
READ(15,*) Nitem
DO i=1, Nitem
   READ(15,*) j, k, m, n, v1, v2, v3, v2BO, vconj
   IF (j == 0 .OR. n==0) THEN
      DO n1 = 1, Natom
         DO n2 = 1, Natom
            torsn(n1,k,m,n2)%V1 = v1
            torsn(n1,k,m,n2)%V2 = v2
            torsn(n1,k,m,n2)%V3 = v3
            torsn(n1,k,m,n2)%V2_BO = v2BO
            torsn(n1,k,m,n2)%vconj = vconj
            IF (k/=m) THEN
               torsn(n2,m,k,n1)%V1    = torsn(n1,k,m,n2)%V1
               torsn(n2,m,k,n1)%V2    = torsn(n1,k,m,n2)%V2
               torsn(n2,m,k,n1)%V3    = torsn(n1,k,m,n2)%V3
               torsn(n2,m,k,n1)%V2_BO = torsn(n1,k,m,n2)%V2_BO
               torsn(n2,m,k,n1)%vconj = torsn(n1,k,m,n2)%vconj
            END IF
         END DO
      END DO
   ELSE
      torsn(j,k,m,n)%V1 = v1
      torsn(j,k,m,n)%V2 = v2
      torsn(j,k,m,n)%V3 = v3
      torsn(j,k,m,n)%V2_BO = v2BO
      torsn(j,k,m,n)%vconj = vconj
      IF (k/=m .OR. j/=n) THEN
         torsn(n,m,k,j)%V1    = torsn(j,k,m,n)%V1
         torsn(n,m,k,j)%V2    = torsn(j,k,m,n)%V2
         torsn(n,m,k,j)%V3    = torsn(j,k,m,n)%V3
         torsn(n,m,k,j)%V2_BO = torsn(j,k,m,n)%V2_BO
         torsn(n,m,k,j)%vconj = torsn(j,k,m,n)%vconj
      END IF
   END IF
END DO
!
! H-Bond potential data parsing
DO i=1, Natom
   hbond(i)%Dehb(:,:) = 0.0D0  ! Criterion for the H-bond calculation
END DO
READ(15,*) Nitem
commn%Npair_HBND = 0
DO i=1, Nitem
   READ(15,*) m, j, n, hbond(j)%Rhb(m,n), hbond(j)%Dehb(m,n), &
        & hbond(j)%vhb1(m,n), hbond(j)%vhb2(m,n)
   commn%Npair_HBND = commn%Npair_HBND + 1
   IF (commn%Npair_HBND > HBOND_max) STOP "===  HBOND max crashed ==="
   commn%HBND_pair(commn%Npair_HBND) = m
   IF (m /= n) THEN
      hbond(j)%Rhb (n,m) = hbond(j)%Rhb (m,n)
      hbond(j)%Dehb(n,m) = hbond(j)%Dehb(m,n)
      hbond(j)%vhb1(n,m) = hbond(j)%vhb1(m,n)
      hbond(j)%vhb2(n,m) = hbond(j)%vhb2(m,n)
      commn%Npair_HBND = commn%Npair_HBND + 1
      IF (commn%Npair_HBND > HBOND_max) STOP "===  HBOND max crashed ==="
      commn%HBND_pair(commn%Npair_HBND) = n
   END IF
END DO
!
param(1)%qmin =  0.0D0; param(1)%qmax =  1.0D0 ! H
param(2)%qmin = -2.0D0; param(2)%qmax =  0.0D0 ! O
param(3)%qmin =  0.0D0; param(3)%qmax =  3.0D0 ! Al
!
RETURN
END SUBROUTINE POTENTIAL_SETUP
