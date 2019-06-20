! function parse(args) result(flags)
!     character(len = 32), dimension(:), intent(in) :: args(:)
!     character(len = 32) :: flags
!     integer :: i

!     do i=1, len(args)
!         if (args(i)(1:1) .eq. "-") then
!             if (args(i)(1:2) .eq. "v") then 
!                 print *, version
!             else if (args(i)(1:2) .eq. "h") then
!                 call print_help()
!             else 
!                 print *, "not a valid flag; run with the ""-h"" flag for options."
!                 stop
!             end if
!         else 
!             print *, "unrecognized input, run with the ""-h"" flag for options."
!         end if
!     end do
    

! end function parse    



program test
! command line argument array
character(len = 32), dimension(:), allocatable :: args(:)
character(len = 32), dimension(:), allocatable :: trimmed(:)
!flag outputs
character(len = 32), parameter :: version = "v1"

integer, dimension(:), allocatable :: blah(:)






integer :: i, num_args
num_args = command_argument_count()
allocate(args(num_args))
allocate(trimmed(num_args))
allocate(blah(0:4))
blah(0) = 1
blah(1) = 1
blah(2) = 2
blah(3) = 3
blah(4) = 4

blah(5) = 4
print *, blah

end program test
! do i = 1, num_args  !populates array of arguments with command line arguments
!     call get_command_argument(i, args(i)) 
! !    write(*, fmt="(1x,a,i0)", advance="no") trim(args(i))
! end do
! do i = 1, num_args
!     trimmed(i) = trim(args(i))
! end do
! !print *, trimmed(1)(1:2)
! call parse(trimmed)

! contains

! subroutine print_help
!      print "(a)", "Usage: ./test [-v -h]" 
!  end subroutine print_help

! subroutine parse(args)
!     character(len = 32), dimension(:), intent(in) :: args(:)   
!     if (len(args) .ne. 0) then
!         do i=1, len(args)
!             if (args(i)(1:1) .eq. "-") then
!                 if (args(i)(2:2) .eq. "v") then 
!                     print *, version
!                     stop
!                 else if (args(i)(2:2) .eq. "h") then
!                     call print_help()
!                     stop
!                 else 
!                     print *, "not a valid flag; run with the ""-h"" flag for options."
!                     stop
!                 end if
!             else 
!                 print *, "unrecognized input, run with the ""-h"" flag for options."
!             end if
!         end do
!     else 
!         print *, "No command line arguments were specified, executing default"
!         stop
!     end if
! end subroutine parse
! end program test