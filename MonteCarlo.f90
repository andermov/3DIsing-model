module Monte_Carlo

public :: Metropolis, neighbours, initialstate
integer, public :: x,y,x2,y2,z,z2,sum1,sum2,sum3,sum4

contains

!This subroutine sets the initial condition of the alloy
!when it is perfectly ordered.
subroutine initialstate(M)

integer, dimension(:,:,:), intent(inout) :: M
integer :: k,l,j,s
real :: n

n=acos(-1.0)
do j=1,20
	do k=1,20
		do l=1,20
			s=k+j
			if (cos(s*n)<0) then
				M(j,k,l)=1
			else
				M(j,k,l)=-1
			endif		
		enddo
	enddo
enddo

end subroutine initialstate

!NB: The subroutine initialstete2 is not used in the execution of the program.
subroutine initialstate2(M)

integer, dimension(:,:,:), intent(inout) :: M
integer :: k,l,j,s
real :: n

call random_seed()
do j=1,20
	do k=1,20
		do l=1,20
			call random_number(n)
			if (n>0.5) then
				M(j,k,l)=1
			else
				M(j,k,l)=-1
			endif		
		enddo
	enddo
enddo

end subroutine initialstate2

subroutine Metropolis(M,p)

integer, dimension(:,:,:), intent(inout) :: M
real, dimension(:), intent (in) :: p
integer :: j,suma,s
real :: n

call random_seed()
do j=1,8000
! Setting a random coordinate in the alloy.
	call random_number(n)
	x=aint(20.0*n+1.0)
	call random_number(n)
	y=aint(20.0*n+1.0)
	call random_number(n)
	z=aint(20.0*n+1.0)
	call random_number(n)
! Interchange in between the previous coordinate and a random neighbour.
	if (n<0.5) then 
		z2=z
		if (n<0.25) then
			x2=x
			if (n<0.125) then
				if (y==20) then
					y2=1
				else
					y2=y+1
				endif
			else
				if (y==1) then
					y2=20
				else
					y2=y-1
				endif
			endif
		else if (n>=0.375) then
			y2=y
			if (x==20) then
				x2=1
			else
				x2=x+1
			endif
		else
			y2=y
			if (x==1) then
				x2=20
			else
				x2=x-1
			endif
		endif
	else
		if (z==20) then
			z2=1
		else
			z2=z+1
		endif
		if (n<0.75) then
			x2=x
			if (n<0.625) then
				if (y==20) then
					y2=1
				else
					y2=y+1
				endif
			else
				if (y==1) then
					y2=20
				else
					y2=y-1
				endif
			endif
		else if (n>=0.875) then
			y2=y
			if (x==20) then
				x2=1
			else
				x2=x+1
			endif
		else
			y2=y
			if (x==1) then
				x2=20
			else
				x2=x-1
			endif
		endif
	endif
	if (M(x,y,z)/=M(x2,y2,z2)) then
!Change in the configuration.
		if (M(x,y,z)==1) then
			call neighbours(x,y,z,suma1,M)
			call neighbours(x2,y2,z2,suma2,M)
			M(x,y,z)=-1
			M(x2,y2,z2)=1
			call neighbours(x,y,z,suma3,M)
			call neighbours(x2,y2,z2,suma4,M)
			suma=-suma1-suma2+suma3+suma4
!Acceptance or rejection of the new configuration.
			if (suma>0) then
				call random_number(n)
				if (n>p(suma)) then
					M(x,y,z)=1
					M(x2,y2,z2)=-1
				endif
			endif
		else
!Change in the configuration.
			call neighbours(x,y,z,suma1,M)
			call neighbours(x2,y2,z2,suma2,M)
			M(x,y,z)=1
			M(x2,y2,z2)=-1
			call neighbours(x,y,z,suma3,M)
			call neighbours(x2,y2,z2,suma4,M)
			suma=-suma1-suma2+suma3+suma4
!Acceptance or rejection of the new configuration.
			if (suma>0) then
				call random_number(n)
				if (n>p(suma)) then
					M(x,y,z)=-1
					M(x2,y2,z2)=1
				endif
			endif
		endif
	endif
enddo

end subroutine Metropolis

!This subroutine counts the number and type of the neighbours.
!Computing the change in the energy depending on the change in the alloy.
subroutine neighbours(x,y,z,intneighbours,M)

integer, dimension(:,:,:), intent(in) :: M
integer, intent(in) :: x,y,z
integer, intent(out) :: intneighbours
integer :: norte, sur, este, oeste, eeste,eoeste,esur,enorte,nneighbours

if (x==1) then
	oeste=M(20,y,z)
else
	oeste=M(x-1,y,z)
endif
if (x==20) then
	este=M(1,y,z)
else
	este=M(x+1,y,z)
endif
if (y==1) then
	sur=M(x,20,z)
else
	sur=M(x,y-1,z)
endif
if (y==20) then
	norte=M(x,1,z)
else
	norte=M(x,y+1,z)
endif
if (z==20) then
	if (x==1) then
		eoeste=M(20,y,1)
	else
		eoeste=M(x-1,y,1)
	endif
	if (x==20) then
		eeste=M(1,y,1)
	else
		eeste=M(x+1,y,1)
	endif
	if (y==1) then
		esur=M(x,20,1)
	else
		esur=M(x,y-1,1)
	endif
	if (y==20) then
		enorte=M(x,1,1)
	else
		enorte=M(x,y+1,1)
	endif
else
	if (x==1) then
		eoeste=M(20,y,z+1)
	else
		eoeste=M(x-1,y,z+1)
	endif
	if (x==20) then
		eeste=M(1,y,z+1)
	else
		eeste=M(x+1,y,z+1)
	endif	
	if (y==1) then
		esur=M(x,20,z+1)
	else
		esur=M(x,y-1,z+1)
	endif
	if (y==20) then
		enorte=M(x,1,z+1)
	else
		enorte=M(x,y+1,z+1)
	endif
endif
nneighbours=(norte+sur+este+oeste+enorte+esur+eoeste+eeste)*M(x,y,z)
intneighbours=nneighbours

end subroutine neighbours

end module Monte_Carlo
