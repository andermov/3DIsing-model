module mc_types

integer, parameter, public :: int4   = selected_int_kind(4)  
integer, parameter, public :: int8   = selected_int_kind(8)  
integer, parameter, public :: int10  = selected_int_kind(10)   

integer, parameter, public :: byte  = int2
integer, parameter, public :: short = int4
integer, parameter, public :: int   = int8
integer, parameter, public :: long  = int10


! Real types

integer, parameter, public :: single = selected_real_kind(6) 
integer, parameter, public :: double = selected_real_kind(14)

! Other useful names for types

integer, parameter, public :: sencillo = single
integer, parameter, public :: doble = double
!
integer, parameter, public :: sp = single
integer, parameter, public :: dp = double

end module mcf_tipos
