program Ising

use mc_types
use Monte_Carlo
use magnitudes

integer, dimension(20,20,20) :: M
real :: T
real (kind=dp) :: Em, Eacum, E2acum, E2m, Hv, E
real, dimension(32) :: p
integer :: d, nequil, i

open (unit=8,status="new", file="100rep20L", action="write")

nequil=100

M=0
call estadoinicial(M)

! Starting from T=0 the temperature increases.
do d=1,400

T=0.05*d
!Boltzmann factors for each temperature.
p(2)=exp(-4.0/T)
p(4)=exp(-8.0/T)
p(6)=exp(-12.0/T)
p(8)=exp(-16.0/T)
p(10)=exp(-20.0/T)
p(12)=exp(-24.0/T)
p(14)=exp(-28.0/T)
p(16)=exp(-32.0/T)
p(18)=exp(-36.0/T)
p(20)=exp(-40.0/T)
p(22)=exp(-44.0/T)
p(24)=exp(-48.0/T)
p(26)=exp(-52.0/T)
p(28)=exp(-56.0/T)
p(30)=exp(-60.0/T)
p(32)=exp(-64.0/T)

!Waiting time until the alloy reaches equilibrium
!after the temperature change.
do i=1,nequil
	call Metropolis(M,p)
enddo

Eacum=0
E2acum=0

!Extracting values for magnitudes of physical significance.
do i=1,100
	call Metropolis(M,p)
	E=0
	call Energia(M,E)
	call datos(Eacum,E2acum,E)
enddo

Em=Eacum/100.0
E2m=E2acum/100.0
Hv=(E2m-Em*Em)/(T**2)
print *, T, Em, Hv
write (unit=8, fmt=*) T, Em, Hv

enddo

end program Ising
