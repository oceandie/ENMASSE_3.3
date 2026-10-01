module types
use par_kind, only: wp        
implicit none
private
public dp, hp, ivector, dvector, zvector

integer, parameter :: dp=wp, &          ! double precision
                      hp=wp             ! high precision

type ivector                       ! allocatable integer vector
   integer, pointer :: vec(:) => null()
end type

type dvector                       ! allocatable real double precision vector
   real(dp), pointer :: vec(:) => null()
end type

type zvector                       ! allocatable complex double precision vector
   complex(dp), pointer :: vec(:) => null()
end type

end module

