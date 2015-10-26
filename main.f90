PROGRAM REAFFOX_AlOH !**********************************************************
!
! REActive Force-Field simulation code for OXidation
!
! ReaxFF library has been given by Dr. Michael Russo at Penn. State Univ.
!
USE DATAFMT
IMPLICIT NONE
! 
TYPE (AMNT):: Nset
TYPE (TCKT):: tag
TYPE (GNRL):: commn
TYPE (ATMC):: param(Natom)
TYPE (ANGL):: angle(Natom,Natom,Natom)
TYPE (HBND):: hbond(Natom)
TYPE (TRSN):: torsn(Natom,Natom,Natom,Natom)
TYPE (STAT):: sys
TYPE (CORR):: pair
REAL   (DP):: dt
INTEGER(GI):: Nloop, OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM, TID, NTD
TYPE(PTCL), POINTER:: q(:)
TYPE(SORT), POINTER:: cell(:)
REAL(SP)   :: secnds, time0, time1, time2
INTEGER(GI):: Nchunk, count_a, count_b, count_rate, count_max
time0 = 0
time1 = secnds(time0)
nset%t_nonb = 0.0D0
nset%t_corr = 0.0D0
nset%t_lpair= 0.0D0
nset%t_over = 0.0D0
nset%t_vac = 0.0D0 
nset%t_tor=0.0D0
nset%t_hbnd = 0.0D0
nset%t_eem1 = 0.0D0
nset%t_eem2 = 0.0D0
nset%t_vv   = 0.0D0
!
!$OMP PARALLEL PRIVATE(TID, NTD)
TID = OMP_GET_THREAD_NUM()
NTD = OMP_GET_NUM_THREADS()
Nset%NTD = NTD
IF (TID == 0) PRINT *, "# Number of available threads = ", NTD
!$OMP END PARALLEL
!
!CALL RANDOM_SEED
CALL POTENTIAL_SETUP(commn, param, angle, hbond, torsn, sys)
CALL READ_INPUT(Nset, tag, sys, q, cell, pair, dt)
CALL CELL_INIT(Nset, cell)
CALL CELL_SORT(Nset, cell, q)
IF (.NOT. tag%restart) &
     & CALL VINIT(Nset, commn, param, angle, hbond, torsn, sys, q, cell, dt)
Nloop = 1
OPEN(UNIT=15, FILE="temporal_status.dat")
WRITE(15,80)
80 FORMAT( "# time(fs) - T(K) - E_total - E_reaxff (eV) - pressure (eV/A3) ")
OPEN(UNIT=17, FILE="e_comp_reaxff.dat")
WRITE(17,90)
90 FORMAT( "# time(fs)  E_bond(2)  E_over(3)  E_under(4) E_lp(5)    ",&
        & "E_val(6)   E_pen(7)   E_H(8)     E_tor(9)   E_conj(10) ", &
        & "E_vdw(11)  E_coul(12) E_charge(13) in kcal/mol")
!
IF (Nset%RUN_type == 0) THEN
   DO WHILE (Nloop <= Nset%freq_max)
      CALL TAG_SET(Nset, tag, Nloop)
      CALL SYSTEM_CLOCK(count_a, count_rate, count_max)
      IF (tag%thermostat .AND. tag%thermo_stoch) THEN
         CALL VV_NVT1_STOCH(Nset, param, q, dt)
      ELSE; CALL VV_NVE1(Nset, param, q, dt); 
      END IF
      CALL SYSTEM_CLOCK(count_b, count_rate, count_max)
      nset%t_vv = nset%t_vv + DBLE(count_b - count_a)/DBLE(count_rate)
      IF (tag%cell_sort)  CALL CELL_SORT(Nset, cell, q)
      IF (tag%new_charge) THEN
         CALL REMOVE_RIGID_MOTION(Nset, q)
         CALL EEM_SHIELDING(Nset, commn, param, sys, cell, q)
      END IF
      CALL FORCE(Nset, commn, param, angle, hbond, torsn, sys, pair, cell, q)
      CALL SYSTEM_CLOCK(count_a, count_rate, count_max)
      IF (tag%thermostat) THEN; 
         IF (tag%thermo_stoch) THEN; CALL VV_NVT2_STOCH(Nset, param, sys, q, dt)
         ELSE; CALL VV_NVT2_BEREND(Nset, param, sys, q, dt)
         END IF
      ELSE; CALL VV_NVE2(Nset, param, sys, q, dt);          
      END IF
      CALL SYSTEM_CLOCK(count_b, count_rate, count_max)
      nset%t_vv = nset%t_vv + DBLE(count_b - count_a)/DBLE(count_rate)
      !
      IF (tag%snapshot)   CALL PRINT_SNAPSHOT(Nset, Nloop, q)   
      IF (tag%print_corr) CALL PRINT_RDF(Nset, Nloop, pair)
      IF (tag%status)     CALL PRINT_STAT(Nset, sys, Nloop, cell, q, dt)
      IF (tag%diff_new)   CALL DIFF_INIT(Nset, q)
      Nloop = Nloop + 1
   END DO
ELSE
   DO WHILE (Nloop <= Nset%freq_max)
      CALL TAG_SET(Nset, tag, Nloop)
      CALL VV_NVE1_REFLECT(Nset, param, sys, q, dt); 
      IF (tag%cell_sort)  CALL CELL_SORT(Nset, cell, q)
      IF (tag%new_charge) THEN
         CALL OXIDE_CHECK(Nset, tag, q)
         IF (tag%new_o2) THEN
            CALL NEW_O2(Nset, sys, q)
            CALL CELL_SORT(Nset, cell, q)
         END IF
         CALL REMOVE_RIGID_MOTION(Nset, q)
         CALL EEM_SHIELDING(Nset, commn, param, sys, cell, q)
      END IF
      CALL FORCE(Nset, commn, param, angle, hbond, torsn, sys, pair, cell, q)
      IF (tag%thermostat) THEN; CALL VV_NVT2_BEREND(Nset, param, sys, q, dt); 
      ELSE; CALL VV_NVE2(Nset, param, sys, q, dt); 
      END IF
      IF (tag%snapshot)   CALL PRINT_SNAPSHOT(Nset, Nloop, q)
      IF (tag%print_corr) CALL PRINT_RDF(Nset, Nloop, pair)
      IF (tag%status)     CALL PRINT_STAT(Nset, sys, Nloop, cell, q, dt)
      IF (tag%diff_new)   CALL DIFF_INIT(Nset, q)
      Nloop = Nloop + 1
      IF (Nset%al_zmax*2.D0 > sys%box(3)) THEN
         PRINT *, "# Oxygen insertion limit crashed"
         GOTO 99
      END IF
   END DO
END IF
!
99 CALL EVACUATE_MEMORY(q, cell)
CLOSE(15)
CLOSE(17)
time2 = secnds(time1)
PRINT '("# Wall time =", ES11.4, " seconds with", I8, " loops")', time2, Nloop-1
PRINT *, "# nonbonded, correction, lone-pair, over/under, vac angle, torsion, hbnd, eem1, eem2, verlet"
PRINT '(10(ES11.4,1X) "sec")', nset%t_nonb, nset%t_corr, nset%t_lpair, &
     & nset%t_over, nset%t_vac, nset%t_tor, nset%t_hbnd, nset%t_eem1, &
     & nset%t_eem2, nset%t_vv
STOP
!
CONTAINS
!
!###############################################################################
SUBROUTINE READ_INPUT(Nset, tag, sys, q, cell, pair, dt)
USE DATAFMT
IMPLICIT NONE
TYPE(AMNT):: Nset
TYPE(TCKT):: tag
TYPE(STAT):: sys
TYPE(CORR):: pair
REAL(DP)  :: dt
TYPE(PTCL), POINTER:: q(:)
TYPE(SORT), POINTER:: cell(:)
!
INTEGER(GI):: openstatus, istatus, i, id, Ncell(3)
REAL   (DP):: qsum, time_max, time_snap, time_rdf, time_stat, time_sort, &
     &        Lcell(3)
CHARACTER(256):: dummy, FILENAME
!
! Initialize tag variables
tag%regular_run = .FALSE.
tag%new_o2      = .FALSE.
tag%new_charge  = .FALSE.
tag%snapshot    = .FALSE.
tag%thermostat  = .FALSE.
tag%restart     = .FALSE.
tag%status      = .FALSE.
tag%print_corr  = .FALSE.
!
! Read simulation parameters **************************************************
OPEN(UNIT=10, FILE='config.aloh', STATUS = "OLD", ACTION="READ", &
     & FORM = "FORMATTED", IOSTAT = openstatus)
IF (openstatus >0) STOP "#### cannot  open config.aloh ###"
READ(10,*) dummy
READ(10,*) dummy
READ(10,*) dummy
CALL ANY_2_UPPER(dummy)
IF (dummy == 'REGULAR') THEN
   Nset%RUN_type = 0
ELSE IF (dummy == 'OXIDATION') THEN
   Nset%RUN_type = 1
ELSE IF (dummy == 'HYDROXIDATION') THEN
   Nset%RUN_type = 2
ELSE
   STOP "=== Simulation type is unknown ==="
END IF
READ(10,*) dummy
READ(10,*) time_max, time_snap, time_rdf, time_stat, time_sort, dt
Nset%freq_max  = NINT(time_max /dt)
Nset%freq_snap = NINT(time_snap/dt)
Nset%freq_rdf  = NINT(time_rdf /dt)
Nset%freq_sort = NINT(time_sort/dt)
Nset%freq_stat = NINT(time_stat/dt)
dt = dt/ TFM
READ(10,*) dummy
READ(10,*) sys%box(:)
Nset%box(:) = sys%box(:)
Nset%V = Nset%box(1)*Nset%box(2)*Nset%box(3)
READ(10,*) dummy
READ(10,*) dummy, sys%T_given, Nset%alpha
sys%T_given = sys%T_given*kB ! K -> eV
CALL ANY_2_UPPER(dummy)
IF (dummy == 'YES') THEN
   tag%thermostat = .TRUE.
ELSE IF (dummy == 'NO') THEN
   tag%thermostat = .FALSE.
ELSE
   STOP "=== Temperature control is unknown ==="
END IF
READ(10,*) dummy
READ(10,*) dummy
CALL ANY_2_UPPER(dummy)
IF (dummy == 'YES') THEN
   tag%restart = .TRUE.
ELSE IF (dummy == 'NO') THEN
   tag%restart = .FALSE.
ELSE
   STOP "=== Restart option is unknown ==="
END IF
sys%Rcut2 = sys%Rcut**2
IF (sys%Rcut*2.D0 > sys%box(1) .OR. sys%Rcut*2.D0 > sys%box(2) .OR. &
& sys%Rcut*2.D0 > sys%box(3)) STOP "=== Rcut is too big ==="
READ(10,*) dummy
READ(10,*) Nset%freq_new   ! New charge frequency
READ(10,*) dummy
READ(10,*) Nset%q_crit     ! New oxygen criterion
READ(10,*) dummy 
READ(10,*) dummy
CALL ANY_2_UPPER(dummy)
IF (dummy == 'YES') THEN
   Nset%tag_msd = .TRUE.
ELSE IF (dummy == 'NO') THEN
   Nset%tag_msd = .FALSE.
ELSE
   STOP "=== MSD print option is unknown ==="
END IF
READ(10,*) dummy
READ(10,*) (Nset%e_field(i),i=1,3) ! E-field
Nset%e_field(:) = Nset%e_field(:)*0.01D0 ! MV/cm => eV/A/e, e=el.charge
CLOSE(10)
IF (Nset%freq_new < 0) THEN
   tag%thermo_stoch = .TRUE.
ELSE
   tag%thermo_stoch = .FALSE.
END IF
!
! Define cell (linked list set) ***********************************************
Ncell(:) = INT(Nset%box(:)/sys%Rcut)
Lcell(:) = Nset%box(:)/DBLE(Ncell(:))
IF (Ncell(1) < 3 .AND. Ncell(2) < 3 .AND. Ncell(3) < 3) THEN
   PRINT *, "==== Error - unitcell is too small or cut-off radius is too large"
   STOP
END IF
WRITE(*,50) Ncell(:), Lcell(:)
50 FORMAT("# === ",I3,"x",I3,"x",I3," cells with size of", 3(F5.1,1X), "===")
Nset%Ncell(:)  = Ncell(:)
Nset%Ncell_all = Ncell(1)*Ncell(2)*Ncell(3)
Nset%Lcell(:)  = Lcell(:)
!
ALLOCATE(cell(Nset%Ncell_all), STAT=istatus)
IF(istatus /=0) STOP "=== Cell allocation error ==="
!
! Read particle data **********************************************************
IF (tag%restart) THEN
   FILENAME = 'restart.xyz'
ELSE
   FILENAME = 'input.xyz'
END IF
OPEN(UNIT=11, FILE=FILENAME, STATUS = "OLD", ACTION="READ", &
     & FORM = "FORMATTED", IOSTAT = openstatus)
IF (openstatus >0) THEN
   PRINT *, '("=== cannot open ", A, "===")', FILENAME; STOP
END IF
!
! Atomic data configuration
! H=1, O=2, Al=3
Nset%Npt(:) = 0
qsum = 0.0D0
Nset%al_zmax = 0.0D0
READ(11,*) Nset%Npt_all
ALLOCATE(q(Nset%Npt_all), STAT=istatus)
IF(istatus /=0) STOP "=== Particle allocation error ==="
READ(11,*) dummy
DO i=1, Nset%Npt_all
   IF (tag%restart) THEN
      READ(11,*) dummy, q(i)%xx(:), q(i)%xv(:), q(i)%ff(:), q(i)%q
      q(i)%xr(:) = q(i)%xx(:)
   ELSE
      READ(11,*) dummy, q(i)%xx(:), q(i)%q
      q(i)%xr(:) = q(i)%xx(:)
   END IF
   q(i)%xx(:) = q(i)%xx(:) - sys%box(:)*DNINT(q(i)%xx(:)/sys%box(:))
   SELECT CASE (dummy)
   CASE ("H")
      id = ID_H
   CASE ("O")
      id = ID_O
   CASE ("Al")
      id = ID_Al
      Nset%al_zmax = MAX(Nset%al_zmax, q(i)%xx(3))
   CASE DEFAULT
      STOP "=== Particle type error ==="
   END SELECT
   Nset%Npt(id) = Nset%Npt(id) + 1
   q(i)%id = id
   qsum = qsum + q(i)%q
END DO
CLOSE(11)
id = 0
DO i=1,Natom
   id = id + Nset%Npt(i)
END DO
IF (id /= Nset%Npt_all) THEN
   DEALLOCATE(q, STAT = istatus)
   STOP "=== Number of particles mismatch ==="
END IF
!
! Charge neutralization
PRINT *, "Charge sum = ", qsum
qsum = qsum/DBLE(Nset%Npt_all)
DO i=1, Nset%Npt_all
   q(i)%q = q(i)%q - qsum
END DO
!
! Correlation data setup ******************************************************
pair%dr = sys%Rcut / DBLE(Nrdf)
pair%V  = Nset%V
pair%rdf(:,:,:) = 0
pair%Npt_old    = Nset%Npt_all
!
RETURN 
END SUBROUTINE READ_INPUT
!
!##############################################################################
SUBROUTINE EVACUATE_MEMORY(q, cell)
USE DATAFMT
IMPLICIT NONE
!
TYPE(PTCL), POINTER:: q(:)
TYPE(SORT), POINTER:: cell(:)
INTEGER(GI):: istatus
DEALLOCATE(q, cell, STAT = istatus)
IF(istatus /=0) STOP "=== Particle/cell deallocation error ==="
!
RETURN
END SUBROUTINE EVACUATE_MEMORY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END PROGRAM REAFFOX_AlOH
!
!http://coding.derkeiler.com/Archive/Fortran/comp.lang.fortran/2005-03/0761.html
!###############################################################################
SUBROUTINE ANY_2_UPPER(txt_string)
USE DATAFMT, ONLY: GI
IMPLICIT NONE
CHARACTER(LEN=*):: txt_string
INTEGER(GI):: i, nlen, id
nlen = LEN(txt_string)
DO i=1, nlen
   id = ichar(txt_string(i:i))
   IF (id >= 97 .AND. id < 122) txt_string(i:i) = CHAR(id-32)
END DO
RETURN
END SUBROUTINE ANY_2_UPPER
!
!###############################################################################
SUBROUTINE TAG_SET(Nset, tag, Nloop)
USE DATAFMT, ONLY: GI, DP, AMNT, TCKT, STAT
IMPLICIT NONE
!
TYPE(AMNT) :: Nset
TYPE(TCKT) :: tag
INTEGER(GI):: Nloop
!
IF (Nset%freq_new < 0) THEN
   tag%new_charge = .FALSE.
ELSE
   IF (MOD(Nloop-1, Nset%freq_new) == 0) THEN
      tag%new_charge = .TRUE.
   ELSE 
      tag%new_charge = .FALSE.
   END IF
END IF
IF (MOD(Nloop, Nset%freq_snap) == 0) THEN
   tag%snapshot = .TRUE.
ELSE
   tag%snapshot = .FALSE.
END IF
IF (MOD(Nloop, Nset%freq_rdf) == 0) THEN
   tag%print_corr = .TRUE.
ELSE
   tag%print_corr = .FALSE.
END IF
IF (MOD(Nloop, Nset%freq_stat) == 0) THEN
   tag%status = .TRUE.
ELSE
   tag%status = .FALSE.
END IF
IF (MOD(Nloop, Nset%freq_sort) == 0) THEN
   tag%cell_sort = .TRUE.
ELSE
   tag%cell_sort = .FALSE.
END IF
RETURN
END SUBROUTINE TAG_SET
