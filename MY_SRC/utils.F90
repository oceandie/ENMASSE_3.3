module utils

! Various general utilities.
! Based on a code by John E. Pask, LLNL.

use types, only: dp
use lib_mpp, only: ctl_stop
implicit none
private
public stop_error

contains

subroutine stop_error(msg)
! Aborts the run on all MPI processes via NEMO's ctl_stop.
! The 'STOP' key makes ctl_stop write the message to ocean.output and call
! mppstop( ld_abort = .true. ), i.e. MPI_ABORT, so the other ranks are killed
! instead of hanging on the next collective communication.
character(len=*), intent(in) :: msg ! Message to print in ocean.output
call ctl_stop( 'STOP', msg )
end subroutine

end module

