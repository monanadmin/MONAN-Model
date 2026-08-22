module NoahmpFatalErrorMod
   ! This module defines an error message for Noah-MP that allows for gracefully failing.
   implicit none
   private
   public  :: NoahMP_error_fatal

   ! I would rather call mpas_derived_types from within this module, but I tried so many
   ! things and I still cannot figure out a way of loading it (somehow mpas_log seems to 
   ! be accessible, but not mpas_derived_types). For the time being I will duplicate the
   ! values, fully aware that this is not the correct thing to do, but please, if you
   ! know how to make this neat, let me know...

   integer, parameter :: NOAHMP_LOG_OUT  = 1  !< code for message type "output"
   integer, parameter :: NOAHMP_LOG_WARN = 2  !< code for message type "warning"
   integer, parameter :: NOAHMP_LOG_ERR  = 3  !< code for message type "error"
   integer, parameter :: NOAHMP_LOG_CRIT = 4  !< code for message type "critical error"

   contains

   !---~---
   !   Print error message then gracefully exit.
   !---~---
   subroutine Noahmp_error_fatal(str)

      use mpas_log          , only : mpas_log_write

      ! input arguments:
      character(len=*),intent(in):: str


      call mpas_log_write(' '                     , messageType=NOAHMP_LOG_ERR)
      call mpas_log_write('---~---'               , messageType=NOAHMP_LOG_ERR)
      call mpas_log_write('   Noah-MP FATAL ERROR', messageType=NOAHMP_LOG_ERR)
      call mpas_log_write('---~---'               , messageType=NOAHMP_LOG_ERR)
      call mpas_log_write(trim(str)               , messageType=NOAHMP_LOG_ERR)
      call mpas_log_write('Noah-MP abort'         , messageType=NOAHMP_LOG_CRIT)
   end subroutine Noahmp_error_fatal
   !---~---
end module NoahmpFatalErrorMod
