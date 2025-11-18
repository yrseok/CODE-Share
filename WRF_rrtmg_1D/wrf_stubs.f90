module module_model_constants
  implicit none
  real, parameter :: cp = 1004.0
end module module_model_constants

module module_state_description
  implicit none
  integer, parameter :: FER_MP_HIRES = 14
  integer, parameter :: FER_MP_HIRES_ADVECT = 15
  integer, parameter :: ETAMP_HWRF = 16
end module module_state_description

module module_wrf_error
contains
  subroutine wrf_error_fatal(msg)
    character(len=*), intent(in) :: msg
    print *, 'FATAL (stub):', trim(msg)
    stop 1
  end subroutine wrf_error_fatal
end module module_wrf_error

logical function wrf_dm_on_monitor()
  wrf_dm_on_monitor = .true.
end function wrf_dm_on_monitor

subroutine wrf_dm_bcast_bytes(buf, nbytes)
  integer, intent(inout) :: buf
  integer, intent(in)    :: nbytes
  ! no-op
end subroutine wrf_dm_bcast_bytes

! 일부 구성에서 호출될 수 있는 CAM 가스 읽기 스텁(안 쓰게 설정하면 호출 안 됨)
module MODULE_RA_CAM_SUPPORT
contains
  subroutine read_CAMgases(yr, julian, plev, tlev, vmr_co2, vmr_ch4, vmr_n2o, vmr_o2, vmr_o3)
    integer, intent(in) :: yr
    real,    intent(in) :: julian
    real,    intent(in) :: plev(:,:), tlev(:,:)
    real,    intent(inout) :: vmr_co2(:,:), vmr_ch4(:,:), vmr_n2o(:,:), vmr_o2(:,:), vmr_o3(:,:)
    ! no-op
  end subroutine read_CAMgases

subroutine wrf_debug(level, msg)
  integer, intent(in) :: level
  character(len=*), intent(in) :: msg
  ! no-op
end subroutine wrf_debug

end module MODULE_RA_CAM_SUPPORT

! =======================
! WRF 전역 스텁 (모듈 바깥)
! =======================

subroutine wrf_error_fatal3(file, line, msg)
  character(len=*), intent(in) :: file, msg
  integer,          intent(in) :: line
  print *, 'FATAL3 (stub):', trim(file), ':', line, ':', trim(msg)
  stop 1
end subroutine wrf_error_fatal3

subroutine wrf_debug(level, msg)
  integer,          intent(in) :: level
  character(len=*), intent(in) :: msg
  ! 원하면 메시지 출력:
  ! print *, 'DEBUG(', level, '): ', trim(msg)
end subroutine wrf_debug

subroutine wrf_message(msg)
  character(len=*), intent(in) :: msg
  ! no-op (필요시 print *, 'MSG: ', trim(msg))
end subroutine wrf_message

