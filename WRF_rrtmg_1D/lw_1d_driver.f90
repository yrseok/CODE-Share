program lw_1d_driver
  use module_ra_rrtmg_lw, only: rrtmg_lwinit, RRTMG_LWRAD 
  use parrrtm,   only: ngptlw
  use parkind,   only: rb => kind_rb
  use rrlw_wvn,  only: rwgt
  implicit none

  character(len=*), parameter :: date_str = "mid-win"

  integer, parameter :: nlev = 33!53!33
  integer, parameter :: nz = nlev !50
  integer, parameter :: nw = nz + 1     ! W-면(계면) 개수

  ! ---- 인덱스(WRF 스타일) ----
  integer :: ids,ide,jds,jde,kds,kde
  integer :: ims,ime,jms,jme,kms,kme
  integer :: its,ite,jts,jte,kts,kte

  ! ---- 상수/제어 ----
  real    :: p_top, psfc, r, g, p0, t_sfc, t_top, lapse
  integer :: icloud, cldovrlp, mp_physics, yr
  logical :: warm_rain, is_cammgmp_used
  integer :: idcor, o3input, ghg_input

  ! ---- 3D (x, k, y): **모두 nw로 통일** ----
  real :: p3d(1,nw,1), t3d(1,nw,1), pi3d(1,nw,1), rho3d(1,nw,1)
  real :: dz8w(1,nw,1)
  real :: qv3d(1,nw,1), qc3d(1,nw,1), qr3d(1,nw,1), qi3d(1,nw,1), qs3d(1,nw,1), qg3d(1,nw,1)
  real :: o33d(1,nw,1)             ! 오존 VMR (mol/mol)

  ! ---- 3D W-면(계면) ----
  real :: p8w(1,nw,1), t8w(1,nw,1)

  ! ---- 2D (x, y) ----
  real :: emiss(1,1), tsk(1,1), xland(1,1), xice(1,1), snow(1,1), xlat(1,1)
  real :: glw(1,1), olr(1,1), lwcf(1,1)

  ! ---- 진단(가열율/플럭스) ----
  real :: rthratenlw(1,nw,1), rthratenlwc(1,nw,1)      ! ★ nw로
  real :: lwupflx(1,nw+1,1), lwupflxc(1,nw+1,1)        ! (= nz+2)
  real :: lwdnflx(1,nw+1,1), lwdnflxc(1,nw+1,1)

  ! ---- 유효반경 & 플래그 ----
  real    :: re_cloud(1,nw,1), re_ice(1,nw,1), re_snow(1,nw,1)
  integer :: has_reqc, has_reqi, has_reqs

  ! ---- 보조 ----
  integer :: i, k
  real :: sigma_m, sigma_w, olr0, t_save

  ! === optical depth arrays (RRTMG LW) ===
  !integer, parameter :: ngptlw = 140      ! RRTMG LW g-points
  !real, dimension(nz+1, ngptlw) :: taut
  !real :: tau_mean(nz)
  
  real(kind=rb) :: taut(nz+1, ngptlw)
  real(kind=rb) :: tautout(nz+1, ngptlw)
  real(kind=rb) :: tau_mean(nz)
  integer :: ib

  ! ---- 출력 CSV 저장 ----
  ! === 추가 선언 (IMPLICIT NONE 바로 아래) ===
  integer :: u
  integer :: u2
  integer :: f
  real, parameter :: sigma = 5.670374419e-8   ! Stefan–Boltzmann constant [W/m^2/K^4]
  real :: B_layer(nz)

  real    :: OLR_base, GLW_base
  real    :: tsk_base, OLR_sfc1, OLR_sfc_small
  real    :: dT_sfc_1, dT_sfc_small
  real    :: dOLR_sfc_1K, dOLR_sfc_small, scale_err
  real    :: dT_layer
  real, dimension(nz) :: dOLR_layer   ! 층별 +1K OLR 변화 저장
  real, dimension(nz) :: dp_layer     ! 각 층의 Δp(가독용, 선택)
  real, dimension(nz) :: contrib_layer     ! 각 층의 Δp(가독용, 선택)
  real, dimension(nz) :: contrib_pct
  real, dimension(nz) :: dOLR_form
  real, dimension(nz) :: contrib_cum

  real, dimension(nz) :: cum_pct, cum_base, cum_base100
  real :: layer_flux(nz), tau_k
  real, dimension(nz) :: lwup_rev, lwdn_rev, p_rev
  real :: K_total, K_atm, cum
  real :: OLR_all1K, dOLR_direct, closure_pct, cum_flux, cum_tot
  real, dimension(nz) :: dLWup, dLWdn, dNet, contrib, cum_pct_tot

  real, dimension(nw) :: lwupflx_base, lwdnflx_base
  real, dimension(nz) :: d_tau

  ! --- 선언부 ---
  real, allocatable :: t3d_save(:,:,:)
  real, allocatable :: tsk_save(:,:)
  real :: p_Pa(nlev), qv(nlev), o3(nlev)

! ================= Profile Setting =================

      !real, parameter :: z1(nlev)=(/   &
      !        0.,    1.,    2.,    3.,    4.,    5.,    6.,    7.,    8.,    &
      !        9.,   10.,   11.,   12.,   13.,   14.,   15.,   16.,   17.,    &  
      !        18.,   19.,   20.,   21.,   22.,   23.,   24.,   25.,   30.,   &
      !        35.,   40.,   45.,   50.,   70.,  100./)
!tropical McClatchey

!      real, parameter :: p1(nlev)=(/ &
!        1.013e+03,9.040e+02,8.050e+02,7.150e+02,6.330e+02,5.590e+02,      &
!        4.920e+02,4.320e+02,3.780e+02,3.290e+02,2.860e+02,2.470e+02,      &
!        2.130e+02,1.820e+02,1.560e+02,1.320e+02,1.110e+02,9.370e+01,      & 
!        7.890e+01,6.660e+01,5.650e+01,4.800e+01,4.090e+01,3.500e+01,      &
!        3.000e+01,2.570e+01,1.220e+01,6.000e+00,3.050e+00,1.590e+00,      & 
!        8.540e-01,5.790e-02,3.000e-04/)
!      real, parameter :: t1(nlev)=(/ &  
!        3.000e+02,2.940e+02,2.880e+02,2.840e+02,2.770e+02,2.700e+02,      & 
!        2.640e+02,2.570e+02,2.500e+02,2.440e+02,2.370e+02,2.300e+02,      & 
!        2.240e+02,2.170e+02,2.100e+02,2.040e+02,1.970e+02,1.950e+02,      &
!        1.990e+02,2.030e+02,2.070e+02,2.110e+02,2.150e+02,2.170e+02,      &
!        2.190e+02,2.210e+02,2.320e+02,2.430e+02,2.540e+02,2.650e+02,      &
!        2.700e+02,2.190e+02,2.100e+02/)
!      real, parameter ::  wh1(nlev)=(/ & 
!             1.900e+01,1.300e+01,9.300e+00,4.700e+00,2.200e+00,1.500e+00,     &
!             8.500e-01,4.700e-01,2.500e-01,1.200e-01,5.000e-02,1.700e-02,     & 
!             6.000e-03,1.800e-03,1.000e-03,7.600e-04,6.400e-04,5.600e-04,     &
!             5.000e-04,4.900e-04,4.500e-04,5.100e-04,5.100e-04,5.400e-04,     &
!             6.000e-04,6.700e-04,3.600e-04,1.100e-04,4.300e-05,1.900e-05,     &
!             6.300e-06,1.400e-07,1.000e-09/)
!      real, parameter ::   wo1(nlev)=(/ & 
!              5.600e-05,5.600e-05,5.400e-05,5.100e-05,4.700e-05,4.500e-05,    &
!              4.300e-05,4.100e-05,3.900e-05,3.900e-05,3.900e-05,4.100e-05,    &
!              4.300e-05,4.500e-05,4.500e-05,4.700e-05,4.700e-05,6.900e-05,    & 
!              9.000e-05,1.400e-04,1.900e-04,2.400e-04,2.800e-04,3.200e-04,    &
!              3.400e-04,3.400e-04,2.400e-04,9.200e-05,4.100e-05,1.300e-05,    & 
!              4.300e-06,8.600e-08,4.300e-11/)
!Mid-summer
!      real, parameter :: p1(nlev)=(/ &
!              1.013e+03,9.020e+02,8.020e+02,7.100e+02,6.280e+02,5.540e+02, &
!              4.870e+02,4.260e+02,3.720e+02,3.240e+02,2.810e+02,2.430e+02, &
!!              2.090e+02,1.790e+02,1.530e+02,1.300e+02,1.110e+02,9.500e+01, &
!              8.120e+01,6.950e+01,5.950e+01,5.100e+01,4.370e+01,3.760e+01, & 
!              3.220e+01,2.770e+01,1.320e+01,6.520e+00,3.330e+00,1.760e+00, &
!              9.510e-01,6.710e-02,3.000e-04/)
!       real, parameter :: t1(nlev)=(/ &
!               2.940e+02,2.900e+02,2.850e+02,2.790e+02,2.730e+02,2.670e+02, &
!               2.610e+02,2.550e+02,2.480e+02,2.420e+02,2.350e+02,2.290e+02, & 
!               2.220e+02,2.160e+02,2.160e+02,2.160e+02,2.160e+02,2.160e+02, &
!               2.160e+02,2.170e+02,2.180e+02,2.190e+02,2.200e+02,2.220e+02, &
!               2.230e+02,2.240e+02,2.340e+02,2.450e+02,2.580e+02,2.700e+02, &
!               2.760e+02,2.180e+02,2.100e+02/)
!       real, parameter ::  wh1(nlev)=(/ &
!               1.400e+01,9.300e+00,5.900e+00,3.300e+00,1.900e+00,1.000e+00, &
!               6.100e-01,3.700e-01,2.100e-01,1.200e-01,6.400e-02,2.200e-02, &
!               6.000e-03,1.800e-03,1.000e-03,7.600e-04,6.400e-04,5.600e-04, &
!               5.000e-04,4.900e-04,4.500e-04,5.100e-04,5.100e-04,5.400e-04, &
!!               6.000e-04,6.700e-04,3.600e-04,1.100e-04,4.300e-05,1.900e-05, &
!               1.300e-06,1.400e-07,1.000e-09/)
!       real, parameter ::   wo1(nlev)=(/ &
!               6.000e-05,6.000e-05,6.000e-05,6.200e-05,6.400e-05,6.600e-05, &
!               6.900e-05,7.500e-05,7.900e-05,8.600e-05,9.000e-05,1.100e-04, &
!               1.200e-04,1.500e-04,1.800e-04,1.900e-04,2.100e-04,2.400e-04, &
!               2.800e-04,3.200e-04,3.400e-04,3.600e-04,3.600e-04,3.400e-04, &
!               3.200e-04,3.000e-04,2.000e-04,9.200e-05,4.100e-05,1.300e-05, &
!               4.300e-06,8.600e-08,4.300e-11/)

!Mid-winter
        real, parameter :: p1(nlev)=(/ &
               1.018e+03,8.973e+02,7.897e+02,6.938e+02,6.081e+02,5.313e+02, &
               4.627e+02,4.016e+02,3.473e+02,2.992e+02,2.568e+02,2.199e+02, &
               1.882e+02,1.610e+02,1.378e+02,1.178e+02,1.007e+02,8.610e+01, &
               7.350e+01,6.280e+01,5.370e+01,4.580e+01,3.910e+01,3.340e+01, &
               2.860e+01,2.430e+01,1.110e+01,5.180e+00,2.530e+00,1.290e+00, & 
               6.820e-01,4.670e-02,3.000e-04/)
       real, parameter :: t1(nlev)=(/ &
               2.722e+02,2.687e+02,2.652e+02,2.617e+02,2.557e+02,2.497e+02, &
               2.437e+02,2.377e+02,2.317e+02,2.257e+02,2.197e+02,2.192e+02, &
               2.187e+02,2.182e+02,2.177e+02,2.172e+02,2.167e+02,2.162e+02, &
               2.157e+02,2.152e+02,2.152e+02,2.152e+02,2.152e+02,2.152e+02, &
               2.152e+02,2.152e+02,2.174e+02,2.278e+02,2.432e+02,2.585e+02, & 
               2.657e+02,2.307e+02,2.102e+02/)
       real, parameter ::  wh1(nlev)=(/ &
               3.500e+00,2.500e+00,1.800e+00,1.200e+00,6.600e-01,3.800e-01, &
               2.100e-01,8.500e-02,3.500e-02,1.600e-02,7.500e-03,6.900e-03, &
               6.000e-03,1.800e-03,1.000e-03,7.600e-04,6.400e-04,5.600e-04, &
               5.000e-04,4.900e-04,4.500e-04,5.100e-04,5.100e-04,5.400e-04, &
               6.000e-04,6.700e-04,3.600e-04,1.100e-04,4.300e-05,1.900e-05, &
               6.300e-06,1.400e-07,1.000e-09/)
       real, parameter ::   wo1(nlev)=(/ &
               6.000e-05,5.400e-05,4.900e-05,4.900e-05,4.900e-05,5.800e-05, &
               6.400e-05,7.700e-05,9.000e-05,1.200e-04,1.600e-04,2.100e-04, &
               2.600e-04,3.000e-04,3.200e-04,3.400e-04,3.600e-04,3.900e-04, &
               4.100e-04,4.300e-04,4.500e-04,4.300e-04,4.300e-04,3.900e-04, &
               3.600e-04,3.400e-04,1.900e-04,9.200e-05,4.100e-05,1.300e-05, &
               4.300e-06,8.600e-08,4.300e-11/)


!gk2a
!      real, parameter ::  p1(nlev)=(/  & 
!              1050.00, 1033.54, 1013.29, 989.45, 962.25, 931.98, 898.92, &
!              863.40, 825.75, 786.32, 745.48, 703.58, 660.99, 618.07, 575.15, &
!              532.57, 490.65, 449.66, 409.88, 371.55, 334.86, 300.00, 267.10, &
!              236.27, 207.60, 181.13, 156.88, 134.83, 114.94, 97.15, 81.37, &
!              67.51, 55.44, 45.04, 36.17, 28.69, 22.46, 17.33, 13.16, 9.83, &
!              7.21, 5.18, 3.64, 2.50, 1.66, 1.08, 0.67, 0.40, 0.23, 0.12, 0.06, &
!              0.03, 0.01/)

!202401040000
!      real, parameter ::  t1(nlev)=(/     &
!              298.9799804688, 298.0799865723, 298.3399963379, 296.2599792480, &
!              294.0299987793, 292.0000000000, 291.3800048828, 289.7699890137, &
!              287.5499877930, 285.6099853516, 284.8599853516, 283.2999877930, &
!              280.9599914551, 278.2999877930, 274.8899841309, 272.2399902344, &
!              268.0199890137, 262.8599853516, 257.9499816895, 252.5800018311, &
!              246.3099975586, 240.6199951172, 234.7299957275, 228.9899902344, &
!              222.5399932861, 216.3699951172, 209.1699981689, 202.6199951172, &
!              195.9700012207, 193.3499908447, 193.4899902344, 197.1999969482, &
!              200.3899993896, 206.7699890137, 209.8699951172, 211.4099884033, &
!              214.5499877930, 218.1199951172, 219.5099945068, 217.4499969482, &
!              217.0499877930, 224.4799957275, 229.9599914551, 237.7399902344, &
!              250.9499969482, 261.6399841309, 260.4400024414, 255.3600006104, &
!!              245.4700012207, 239.9299926758, 230.9299926758, 215.8799896240, &
!              183.3699951172/)
              
! 202001100000
!      real, parameter ::  t1(nlev)=(/     &
!              298.6399841309, 297.7500000000, 296.3500061035, 295.6099853516, &
!              293.4299926758, 291.9199829102, 290.8299865723, 289.3099975586, &
!              287.5599975586, 285.0499877930, 283.1999816895, 281.8699951172, &
!              279.3399963379, 277.2999877930, 272.6600036621, 267.8200073242, &
!              263.2299804688, 258.7900085449, 253.7799987793, 249.1799926758, &
!              244.3699951172, 239.1499938965, 235.1999969482, 230.9899902344, &
!              225.5499877930, 218.8799896240, 212.5599975586, 204.7699890137, &
!              198.2399902344, 195.1699981689, 194.8399963379, 195.4400024414, &
!             202.9899902344, 207.9799957275, 213.5899963379, 214.9299926758, &
!              218.3999938965, 222.2999877930, 224.6100006104, 226.7200012207, &
!              227.8300018311, 229.2299957275, 235.2399902344, 249.3699951172, &
!              258.4199829102, 262.5400085449, 262.7699890137, 259.6900024414, &
!              247.6999969482, 239.7999877930, 226.9599914551, 213.5199890137, &
!              184.6799926758/)
!202401040000
!            real, parameter :: wh1(nlev)=(/ &
!              15.2299995422, 15.2299995422, 15.1899995804, 14.9600000381, &
!              14.6700000763, 13.8099994659, 11.6799993515, 10.2199993134, &
!              9.3599996567, 7.2799997330, 2.1199998856, 0.8100000024, &
!              1.9399999380, 3.3599998951, 2.4299998283, 1.7199999094, &
!              0.9300000072, 0.3499999940, 0.2800000012, 0.4899999797, &
!              0.2699999809, 0.0499999970, 0.0299999993, 0.0199999996, &
!              0.0099999998, 0.0000000000, 0.0000000000, 0.0000000000, &
!              0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, &
!              0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, &
!              0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, &
!              0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, &
!              0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, &
!              0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, &
!              0.0000000000/)


! 202001100000
!      real, parameter :: wh1(nlev)=(/ &
!              15.4499998093, 15.4499998093, 15.4499998093, 15.0799999237, 14.5099992752, &
!              12.8800001144, 11.0500001907, 9.4799995422, 6.7500000000, 5.1099996567, &
!             2.1599998474, 1.0199999809, 0.6200000048, 0.3499999940, 0.1999999881, &
!              0.3799999952, 0.1899999976, 0.3599999845, 0.5399999619, 0.5299999714, &
!              0.4599999785, 0.3100000024, 0.0599999987, 0.0000000000, 0.0000000000, &
!              0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, &
!              0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, &
!              0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, &
!              0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, &
!              0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, &
!              0.0000000000, 0.0000000000, 0.0000000000/)
              

!      real, parameter ::  wo1(nlev)=(/ &
!              0.00 , 0.000, 0.00, 0.00, 0.00, 0.000, 0.00, 0.00, 0.00, 0.00, 0.000, 0.00, &
!              0.000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, &
!              0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,&
!!              0.0000000000, 0.0000000000, 0.0000000000, 0.00 , 0.000, 0.00, 0.00, &
!              0.00, 0.000, 0.00, 0.00, 0.00, 0.00, 0.000, 0.00, &
!              0.000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, &
!              0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,&
!              0.0000000000, 0.0000000000, 0.0000000000, 0.00/)

! 단위 변환
p_Pa   = p1 * 100.0       ! hPa -> Pa
qv   = wh1 * 1.0e-3     ! g/kg -> kg/kg
o3   = wo1              ! [mol/mol], 그대로 사용
! --- 기본 경계값 업데이트 ---
p_top = max(minval(p_Pa),100.0)
psfc  = maxval(p_Pa)
t_sfc = t1(1)
t_top = t1(nlev)
! --- 지표 조건 (필수)
icloud     = 0
o3input    = 1
ghg_input  = 0
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ids=1; ide=1; jds=1; jde=1
  ims=1; ime=1; jms=1; jme=1
  its=1; ite=1; jts=1; jte=1
  kds=1; kde=nw         ! ★ 도메인 k-상한 = nw
  kms=1; kme=nw         ! ★ 메모리 k-상한 = nw (더미 요소 포함)
  kts=1; kte=nz         ! ★ 연산은 질량면 nz까지 (p8w는 내부에서 kte+1 접근)


! --- RRTMG 입력 배열 구성 (top→bottom 순서 유지) ---
do k = 1, nz
  p3d(1,k,1)   = p_Pa(k)
  t3d(1,k,1)   = T1(k)
  qv3d(1,k,1)  = qv(k) * 1.0e-3
  o33d(1,k,1)  = o3(k)
  rho3d(1,k,1) = p_Pa(k) / (r * T1(k))
  pi3d(1,k,1)  = (p_Pa(k)/p0)**(r/1004.0)
end do


! --- 상단 더미층 (유지)
p3d(1,nw,1)=p3d(1,nz,1)
t3d(1,nw,1)=t3d(1,nz,1)
rho3d(1,nw,1)=rho3d(1,nz,1)
pi3d(1,nw,1)=pi3d(1,nz,1)
dz8w(1,nw,1)=0.0
qv3d(1,nw,1)=0.0
o33d(1,nw,1)=0.0

re_cloud(1,nw,1)=re_cloud(1,nz,1)
re_ice(1,nw,1)  =re_ice(1,nz,1)
re_snow(1,nw,1) =re_snow(1,nz,1)
print *, 'cloud: ', re_cloud(1,1,1)
print *, 'cloud: ', re_cloud(1,2,1)
! ================= 계면(W-면) =================
do k=1,nw
  sigma_w     = real(k-1)/real(nw-1)            ! 0..1 (지표→꼭대기)
  p8w(1,k,1)  = psfc - (psfc - p_top)*sigma_w    ! ← 방향 수정
!  p8w(1,k,1)  = p_top + (psfc - p_top)*sigma_w  ! top→bottom 순서
end do

! W-면 온도
t8w(1,1,1)=t3d(1,1,1)
do k=2,nz
  t8w(1,k,1)=0.5*(t3d(1,k-1,1)+t3d(1,k,1))
end do
t8w(1,nw,1)=t3d(1,nz,1)

  ! ================= 표면/위치 =================
  emiss(1,1)=1.0;  tsk(1,1)=t_sfc
  xland(1,1)=1.0;  xice(1,1)=0.0;  snow(1,1)=0.0
  xlat(1,1)=0.0

  has_reqc=0; has_reqi=0; has_reqs=0
  rthratenlw=0.0; rthratenlwc=0.0
  lwupflx=0.0; lwupflxc=0.0; lwdnflx=0.0; lwdnflxc=0.0

  ! ================= 초기화 =================
  call rrtmg_lwinit( p_top, .true., &
                     ids,ide,jds,jde,kds,kde, &
                     ims,ime,jms,jme,kms,kme, &
                     its,ite,jts,jte,kts,kte )

  ! ================= Baseline =================
  print *, 'tile i:', its, ite, '  tile j:', jts, jte, '  k:', kts, kte
  print *, 'p8w(1,1,1)=', p8w(1,1,1), '  tsk=', tsk(1,1), '  t3d(1,1,1)=', t3d(1,1,1)
  print *, 'p8w(1,1,1)=', p8w(1,1,1), '  p8w(1,nw,1)=', p8w(1,nw,1)
  print *, 'nw= ', nw

  ! 기대: p8w(1,1,1) ≈ psfc(=1e5 Pa) > p8w(1,nw,1) ≈ p_top(=5e3 Pa)

! 단조 감소 검증 (위로 갈수록 압력 감소)
do k=1,nz
  if ( p8w(1,k,1) <= p8w(1,k+1,1) ) then
    print *, 'ERROR: p8w not strictly decreasing at k=', k, p8w(1,k,1), p8w(1,k+1,1)
    stop 2
  end if
end do

!============================================================
! --- Normalize RRTMG Gaussian weights (only once) ---
rwgt = rwgt !/ sum(rwgt)
! --- Save File CSV ---
open(newunit=u, file="rwgt_" // date_str // ".csv", status='replace', action='write')
write(u,'(A)') 'rwgt'

do k = 1, ngptlw
   write(u,'(F10.6)') &
        rwgt(k)
end do

close(u)

rwgt = rwgt / sum(rwgt)
!============================================================
! (놓치기 쉬움) 단일 칼럼 타일 인덱스
its=1; ite=1
jts=1; jte=1

  yr=2020; glw=0.0; olr=0.0; lwcf=0.0
  call RRTMG_LWRAD( &
    rthratenlw=rthratenlw, rthratenlwc=rthratenlwc, &
    glw=glw, olr=olr, lwcf=lwcf, emiss=emiss, &
    p8w=p8w, p3d=p3d, pi3d=pi3d, dz8w=dz8w, &
    tsk=tsk, t3d=t3d, t8w=t8w, rho3d=rho3d, r=r, g=g, &
    icloud=icloud, warm_rain=warm_rain, cldovrlp=cldovrlp, &
    is_cammgmp_used=is_cammgmp_used, &
    xland=xland, xice=xice, snow=snow, xlat=xlat, &
    qv3d=qv3d, qc3d=qc3d, qr3d=qr3d, qi3d=qi3d, qs3d=qs3d, qg3d=qg3d, &
    re_cloud=re_cloud, re_ice=re_ice, re_snow=re_snow, &
    has_reqc=has_reqc, has_reqi=has_reqi, has_reqs=has_reqs, &
    calc_clean_atm_diag=0, yr=yr, julian=180.0, mp_physics=mp_physics, &
    ids=ids,ide=ide, jds=jds,jde=jde, kds=kds,kde=kde, &
    ims=ims,ime=ime, jms=jms,jme=jme, kms=kms,kme=kme, &
    its=its,ite=ite, jts=jts,jte=jte, kts=kts,kte=kte, &
    lwupflx=lwupflx, lwupflxc=lwupflxc, lwdnflx=lwdnflx, lwdnflxc=lwdnflxc, &
    idcor=idcor, o3input=o3input, o33d=o33d, ghg_input=ghg_input, taut=tautout )

  olr0 = olr(1,1)
  write(*,'(A,F10.3)') 'Baseline OLR [W/m^2] = ', olr0
  write(*,'(A,F10.3)') 'Baseline GLW [W/m^2] = ', glw(1,1) !Ground Downwelling Longwave Rad

  write(*,'(A,F10.3)') 'lwdnflx [W/m^2] = ', lwdnflx(1,4,1) !Ground Downwelling Longwave Rad
  write(*,'(A,F10.3)') 'lwdnflxc [W/m^2] = ', lwdnflxc(1,4,1) !Ground Downwelling Longwave Rad


do k = 1, nw
   lwupflx_base(k) = lwupflx(1,k,1)
   lwdnflx_base(k) = lwdnflx(1,k,1)
!   print *, 'lwupflx_base = ', lwupflx(1,k,1)
!   print *, 'lwdnflx_base = ', lwdnflx(1,k,1)

enddo

do k = 1, nz
    ! 상향 플럭스와 하향 플럭스 차이를 사용하여 Δτ 계산
    d_tau(k) = log( (lwdnflx(1,k+1,1) - lwupflx(1,k,1)) / &
                        (lwdnflx(1,k,1) - lwupflx(1,k+1,1)) )
!    print *, 'lwdnflx_next = ', lwdnflx(1,k+1,1)
!    print *, 'lwupflx_next = ', lwupflx(1,k+1,1)
!    print *, 'lwdnflx_base = ', lwdnflx(1,k,1)
!    print *, 'lwupflx_base = ', lwupflx(1,k,1)
!    print *, 'exp = ', exp(d_tau(k))
!    print *, 'd_tau = ', d_tau(k)
end do
do k = 1, nz
   tau_mean(k) = 0.0
   do ib = 1, ngptlw
      tau_mean(k) = tau_mean(k) + rwgt(ib) * exp(-tautout(k, ib))
   end do
   tau_mean(k) = -log(tau_mean(k))
end do

B_layer = 0.0
! --- 각 층의 Planck 복사량 계산 ---
do k = 1, nz
   B_layer(k) = sigma * t3d(1,k,1)**4
end do

! --- Save File CSV ---
open(newunit=u, file="olr_mclr_c1_" // date_str // ".csv", status='replace', action='write')

! 헤더를 작성 (각 변수의 이름)
write(u,'(A)') 'k,p_mid_Pa,T_mid_K,lwup_Wm2,lwdn_Wm2,lwup_next_Wm2,lwdn_next_Wm2,d_tau,B_layer_Wm2,tau_mean'

! 각 층에 대해 필요한 정보 저장
do k = 1, nz
   write(u,'(I4,1x,F12.3,1x,F9.3,1x,E12.5,1x,E12.5,1x,E12.5,1x,E12.5,1x,E12.5x,E12.5x,E12.5)') & 
        k, p3d(1,k,1), t3d(1,k,1), lwupflx(1,k,1), lwdnflx(1,k,1), &
        lwupflx(1,k+1,1), lwdnflx(1,k+1,1), d_tau(k), B_layer(k), tau_mean(k)
end do
close(u)
print *, '-> wrote olr_kernels_mclr_c1_' // date_str // '.csv'

! g-point tau per layer
open(newunit=u2, file="olr_tau_c1_" // date_str // ".csv", status='replace', action='write')
! ---- Header ----
write(u2,'(A)', advance='no') 'layer'
do f = 1, ngptlw
    write(u2,'(A,I3)', advance='no') ',tau_g', f
end do
write(u2,*)    ! finish header line

! ---- Rows: each layer has 140 tau values horizontally ----
do k = 1, nz
    write(u2,'(I4)', advance='no') k   ! layer index first

    do f = 1, ngptlw
        write(u2,'(",",ES12.6)', advance='no') taut(k, f)
    end do

    write(u2,*)     ! end of the row
end do

close(u2)
print *, '-> wrote ', 'c1_g-point'


! ===== CASE 2 Surface Kernel (+1K) & 선형성(+0.1K) 체크 =====
tsk_base = tsk(1,1)

! --- +1 K ---
tsk(1,1) = tsk_base + 1.0
call RRTMG_LWRAD( &
  rthratenlw=rthratenlw, rthratenlwc=rthratenlwc, &
  glw=glw, olr=olr, lwcf=lwcf, emiss=emiss, &
  p8w=p8w, p3d=p3d, pi3d=pi3d, dz8w=dz8w, &
  tsk=tsk, t3d=t3d, t8w=t8w, rho3d=rho3d, r=r, g=g, &
  icloud=icloud, warm_rain=warm_rain, cldovrlp=cldovrlp, &
  is_cammgmp_used=is_cammgmp_used, &
  xland=xland, xice=xice, snow=snow, xlat=xlat, &
  qv3d=qv3d, qc3d=qc3d, qr3d=qr3d, qi3d=qi3d, qs3d=qs3d, qg3d=qg3d, &
  re_cloud=re_cloud, re_ice=re_ice, re_snow=re_snow, &
  has_reqc=has_reqc, has_reqi=has_reqi, has_reqs=has_reqs, &
  calc_clean_atm_diag=0, yr=yr, julian=180.0, mp_physics=mp_physics, &
  ids=ids,ide=ide, jds=jds,jde=jde, kds=kds,kde=kde, &
  ims=ims,ime=ime, jms=jms,jme=jme, kms=kms,kme=kme, &
  its=its,ite=ite, jts=jts,jte=jte, kts=kts,kte=kte, &
  lwupflx=lwupflx, lwupflxc=lwupflxc, lwdnflx=lwdnflx, lwdnflxc=lwdnflxc, &
  idcor=idcor, o3input=o3input, o33d=o33d, ghg_input=ghg_input, taut=tautout )
dOLR_sfc_1K = olr(1,1) - olr0
print *, 'ΔOLR(surface +1 K) = ', dOLR_sfc_1K


B_layer = 0.0
! --- 각 층의 Planck 복사량 계산 ---
do k = 1, nz
   B_layer(k) = sigma * t3d(1,k,1)**4
end do

! --- Save File CSV ---
open(newunit=u, file="olr_mclr_c2_" // date_str // ".csv", status='replace', action='write')
write(u,'(A)') 'k,p_mid_Pa,T_mid_K,lwup_Wm2,B_layer_Wm2'

do k = 1, nz
   write(u,'(I4,1x,F12.3,1x,F9.3,1x,E12.5,1x,E10.3)') &
        k, p3d(1,k,1), t3d(1,k,1), lwupflx(1,k,1), B_layer(k)
end do

close(u)
print *, '-> wrote olr_kernels_mclr_c2_' // date_str // '.csv'


open(newunit=u2, file="olr_tau_c2_" // date_str // ".csv", status='replace', action='write')
! ---- Header ----
write(u2,'(A)', advance='no') 'layer'
do f = 1, ngptlw
    write(u2,'(A,I3)', advance='no') ',tau_g', f
end do
write(u2,*)    ! finish header line

! ---- Rows: each layer has 140 tau values horizontally ----
do k = 1, nz
    write(u2,'(I4)', advance='no') k   ! layer index first

    do f = 1, ngptlw
        write(u2,'(",",ES12.6)', advance='no') taut(k, f)
    end do

    write(u2,*)     ! end of the row
end do

close(u2)
print *, '-> wrote ', 'c2_g-point'


! only surface +1 K (CASE2) -------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! --- +0.1 K (선형성 확인용) ---
tsk(1,1) = tsk_base + 0.1
call RRTMG_LWRAD( &
  rthratenlw=rthratenlw, rthratenlwc=rthratenlwc, &
  glw=glw, olr=olr, lwcf=lwcf, emiss=emiss, &
  p8w=p8w, p3d=p3d, pi3d=pi3d, dz8w=dz8w, &
  tsk=tsk, t3d=t3d, t8w=t8w, rho3d=rho3d, r=r, g=g, &
  icloud=icloud, warm_rain=warm_rain, cldovrlp=cldovrlp, &
  is_cammgmp_used=is_cammgmp_used, &
  xland=xland, xice=xice, snow=snow, xlat=xlat, &
  qv3d=qv3d, qc3d=qc3d, qr3d=qr3d, qi3d=qi3d, qs3d=qs3d, qg3d=qg3d, &
  re_cloud=re_cloud, re_ice=re_ice, re_snow=re_snow, &
  has_reqc=has_reqc, has_reqi=has_reqi, has_reqs=has_reqs, &
  calc_clean_atm_diag=0, yr=yr, julian=180.0, mp_physics=mp_physics, &
  ids=ids,ide=ide, jds=jds,jde=jde, kds=kds,kde=kde, &
  ims=ims,ime=ime, jms=jms,jme=jme, kms=kms,kme=kme, &
  its=its,ite=ite, jts=jts,jte=jte, kts=kts,kte=kte, &
  lwupflx=lwupflx, lwupflxc=lwupflxc, lwdnflx=lwdnflx, lwdnflxc=lwdnflxc, &
  idcor=idcor, o3input=o3input, o33d=o33d, ghg_input=ghg_input,taut=tautout )
  dOLR_sfc_small = olr(1,1) - olr0

scale_err = 100.0 * abs(dOLR_sfc_1K - 10.0*dOLR_sfc_small) / max(1.0, abs(dOLR_sfc_1K))
tsk(1,1) = tsk_base   ! 표면 온도 원복

OLR_base = OLR(1,1)
!write(*,'(A,F10.3)') 'Baseline TOA OLR [W/m^2] = ', OLR_base

!write(*,'(A,F10.3)')            'dOLR_surface(+1K) [W/m^2] = ', dOLR_sfc_1K
!write(*,'(A,F10.3,A,F6.2)') 'dOLR_surface(+0.1K)*10 = ', 10.0*dOLR_sfc_small, '   scale_err(%)=', scale_err
!print *, 'lwupflx(1,1,1)=', lwupflx(1,1,1)
!print *, 'lwupflx(1,nw,1)=', lwupflx(1,nw,1)
!print *, 'OLR(1,1)=', OLR(1,1)


  ! ================= per layer +1K =================
  do k=1,nz
    t_save      = t3d(1,k,1)
    t3d(1,k,1)  = t_save + 1.0
    t8w(1,1,1)  = t3d(1,1,1)
    do i=2,nz
      t8w(1,i,1)= 0.5*(t3d(1,i-1,1)+t3d(1,i,1))
    end do
    t8w(1,nw,1) = t3d(1,nz,1)

    glw=0.0; olr=0.0; lwcf=0.0
    rthratenlw=0.0; rthratenlwc=0.0
    lwupflx=0.0; lwupflxc=0.0; lwdnflx=0.0; lwdnflxc=0.0

    call RRTMG_LWRAD( &
      rthratenlw=rthratenlw, rthratenlwc=rthratenlwc, &
      glw=glw, olr=olr, lwcf=lwcf, emiss=emiss, &
      p8w=p8w, p3d=p3d, pi3d=pi3d, dz8w=dz8w, &
      tsk=tsk, t3d=t3d, t8w=t8w, rho3d=rho3d, r=r, g=g, &
      icloud=icloud, warm_rain=warm_rain, cldovrlp=cldovrlp, &
      is_cammgmp_used=is_cammgmp_used, &
      xland=xland, xice=xice, snow=snow, xlat=xlat, &
      qv3d=qv3d, qc3d=qc3d, qr3d=qr3d, qi3d=qi3d, qs3d=qs3d, qg3d=qg3d, &
      re_cloud=re_cloud, re_ice=re_ice, re_snow=re_snow, &
      has_reqc=has_reqc, has_reqi=has_reqi, has_reqs=has_reqs, &
      calc_clean_atm_diag=0, yr=yr, julian=180.0, mp_physics=mp_physics, &
      ids=ids,ide=ide, jds=jds,jde=jde, kds=kds,kde=kde, &
      ims=ims,ime=ime, jms=jms,jme=jme, kms=kms,kme=kme, &
      its=its,ite=ite, jts=jts,jte=jte, kts=kts,kte=kte, &
      lwupflx=lwupflx, lwupflxc=lwupflxc, lwdnflx=lwdnflx, lwdnflxc=lwdnflxc, &
      idcor=idcor, o3input=o3input, o33d=o33d, ghg_input=ghg_input,taut=tautout )

!    write(*,'(A,I3,A,F10.3)') 'dOLR for layer ', k, ' = ', olr(1,1)-olr0
    t3d(1,k,1)  = t_save
    dOLR_layer(k) = olr(1,1) - olr0
    dp_layer(k)   = p8w(1,k,1) - p8w(1,k+1,1)

  end do


B_layer = 0.0
! --- 각 층의 Planck 복사량 계산 ---
do k = 1, nz
   B_layer(k) = sigma * t3d(1,k,1)**4
end do

! --- Save File CSV ---
open(newunit=u, file="olr_mclr_c3_" // date_str // ".csv", status='replace', action='write')
write(u,'(A)') 'k,p_mid_Pa,T_mid_K,lwup_Wm2,B_layer_Wm2'

do k = 1, nz
   write(u,'(I4,1x,F12.3,1x,F9.3,1x,E12.5,1x,E10.3)') &
        k, p3d(1,k,1), t3d(1,k,1), lwupflx(1,k,1), B_layer(k)
end do

close(u)
print *, '-> wrote olr_kernels_mclr_c3_' // date_str // '.csv'


open(newunit=u2, file="olr_tau_c3_" // date_str // ".csv", status='replace', action='write')
! ---- Header ----
write(u2,'(A)', advance='no') 'layer'
do f = 1, ngptlw
    write(u2,'(A,I3)', advance='no') ',tau_g', f
end do
write(u2,*)    ! finish header line

! ---- Rows: each layer has 140 tau values horizontally ----
do k = 1, nz
    write(u2,'(I4)', advance='no') k   ! layer index first

    do f = 1, ngptlw
        write(u2,'(",",ES12.6)', advance='no') taut(k, f)
    end do

    write(u2,*)     ! end of the row
end do

close(u2)
print *, '-> wrote ', 'c3_g-point'



!################## kernel
K_atm  = sum(dOLR_layer(1:nz))
K_total = dOLR_sfc_1K + K_atm

!cum = dOLR_sfc_1K
cum = 0.0
do k = 1, nz
  !cum = cum + dOLR_layer(k)
  !cum_pct(k) = 100.0 * cum / max(1.0e-8, K_total)
  cum = cum + dOLR_layer(k)*dp_layer(k)
  cum_pct(k) = 100.0 * cum / sum(dOLR_layer(1:nz) * dp_layer(1:nz))
  !write(*,'(I3,1X,F10.2,1X,F10.4,1X,F8.2)') k, p3d(1,k,1), dOLR_layer(k), cum_pct(k)

end do

write(*,'(A,F7.3)') 'Sum of layer kernels (atm)  [W/m^2/K] (K_atm) = ', K_atm
write(*,'(A,F7.3)') 'Total kernel (atm+surface) [W/m^2/K] (K_total) = ', K_total

! ------------------------------
! 전층+표면 +1K 직접 테스트 (Closure)
! ------------------------------
! --- 배열 복사용 변수 준비
allocate(t3d_save(size(t3d,1), size(t3d,2), size(t3d,3)))
allocate(tsk_save(size(tsk,1), size(tsk,2)))

t3d_save = t3d
tsk_save = tsk

! --- 전층+표면 온도 +1K ---
do k = 1, nz
  t3d(1,k,1) = t3d(1,k,1) + 1.0
end do
tsk(1,1) = tsk(1,1) + 1.0

call RRTMG_LWRAD( &
     rthratenlw=rthratenlw, rthratenlwc=rthratenlwc, &
     glw=glw, olr=olr, lwcf=lwcf, emiss=emiss, &
     p8w=p8w, p3d=p3d, pi3d=pi3d, dz8w=dz8w, &
     tsk=tsk, t3d=t3d, t8w=t8w, rho3d=rho3d, r=r, g=g, &
     icloud=icloud, warm_rain=warm_rain, cldovrlp=cldovrlp, &
     is_cammgmp_used=is_cammgmp_used, &
     xland=xland, xice=xice, snow=snow, xlat=xlat, &
     qv3d=qv3d, qc3d=qc3d, qr3d=qr3d, qi3d=qi3d, qs3d=qs3d, qg3d=qg3d, &
     re_cloud=re_cloud, re_ice=re_ice, re_snow=re_snow, &
     has_reqc=has_reqc, has_reqi=has_reqi, has_reqs=has_reqs, &
     calc_clean_atm_diag=0, yr=yr, julian=180.0, mp_physics=mp_physics, &
     ids=ids,ide=ide, jds=jds,jde=jde, kds=kds,kde=kde, &
     ims=ims,ime=ime, jms=jms,jme=jme, kms=kms,kme=kme, &
     its=its,ite=ite, jts=jts,jte=jte, kts=kts,kte=kte, &
     lwupflx=lwupflx, lwupflxc=lwupflxc, lwdnflx=lwdnflx, lwdnflxc=lwdnflxc, &
     idcor=idcor, o3input=o3input, o33d=o33d, ghg_input=ghg_input,taut=tautout )

OLR_all1K   = OLR(1,1)
dOLR_direct = OLR_all1K - OLR_base

B_layer = 0.0
! --- 각 층의 Planck 복사량 계산 ---
do k = 1, nz
   B_layer(k) = sigma * t3d(1,k,1)**4
end do

! --- Save File CSV ---
open(newunit=u, file="olr_mclr_c4_" // date_str // ".csv", status='replace', action='write')
write(u,'(A)') 'k,p_mid_Pa,T_mid_K,lwup_Wm2,B_layer_Wm2'

do k = 1, nz
   write(u,'(I4,1x,F12.3,1x,F9.3,1x,E12.5,1x,E10.3)') &
        k, p3d(1,k,1), t3d(1,k,1), lwupflx(1,k,1), B_layer(k)
end do

close(u)
print *, '-> wrote olr_kernels_mclr_c4_' // date_str // '.csv'


open(newunit=u2, file="olr_tau_c4_" // date_str // ".csv", status='replace', action='write')
! ---- Header ----
write(u2,'(A)', advance='no') 'layer'
do f = 1, ngptlw
    write(u2,'(A,I3)', advance='no') ',tau_g', f
end do
write(u2,*)    ! finish header line

! ---- Rows: each layer has 140 tau values horizontally ----
do k = 1, nz
    write(u2,'(I4)', advance='no') k   ! layer index first

    do f = 1, ngptlw
        write(u2,'(",",ES12.6)', advance='no') taut(k, f)
    end do

    write(u2,*)     ! end of the row
end do

close(u2)
print *, '-> wrote ', 'c4_g-point'


print *, "τ_mean =", tau_mean(1)
print *, "τ(1,1) =", taut(1,1)
print *, 'DEBUG: taut(1,1)=', taut(1,1)
print *, 'DEBUG: taut(1,10)=', taut(1,10)
print *, 'DEBUG: taut(10,1)=', taut(10,1)
print *, 'shape(taut) = ', size(taut,1), size(taut,2)
print *, 'sum(rwgt) =', sum(rwgt)
print *, 'min(rwgt), max(rwgt) =', minval(rwgt), maxval(rwgt)

!print *, '=============================================='
!print *, ' Layer  k   G-point  τ(k,ib)'
!print *, '=============================================='
!do k = 1, 1
!   do ib = 1, 140   ! ngptlw
!!      write(*,'(I4,1X,I4,1X,E12.5)') k, ib, taut(k,ib)
!   end do
!end do
!print *, '=============================================='

! --- 원복 & 해제 ---
t3d = t3d_save
tsk = tsk_save
deallocate(t3d_save)

contains
! ===== 로그-보간 함수 추가 =====
real function log_interp(xarr, yarr, xval, n)
  implicit none
  integer, intent(in) :: n
  real, intent(in) :: xarr(n), yarr(n), xval
  integer :: i
  if (xval <= xarr(n)) then
    log_interp = yarr(n)
  else if (xval >= xarr(1)) then
    log_interp = yarr(1)
  else
    do i=1,n-1
      if (xval <= xarr(i) .and. xval >= xarr(i+1)) then
        log_interp = yarr(i) + (yarr(i+1)-yarr(i))*(xval-xarr(i))/(xarr(i+1)-xarr(i))
        return
      end if
    end do
  end if
end function log_interp


end program lw_1d_driver

