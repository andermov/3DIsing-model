module magnitudes

use mcf_tipos

public :: Energia, datos

contains

subroutine Energia(M,E)

integer, dimension(:,:,:), intent(in) :: M
real (kind=dp), intent (inout) :: E
integer :: j,k,l,suma,norte,este,enorte,eeste

do j=1,20
	do k=1,20
		do l=1,20
			if (M(j,k,l)==1) then
				if (k==20) then
					norte=M(j,1,l)
				else
					norte=M(j,k+1,l)
				endif
				if (j==20) then
					este=M(1,k,l)
				else
					este=M(j+1,k,l)
				endif
				if (l==20) then
					if (j==20) then
						eeste=M(1,k,1)
					else
						eeste=M(j+1,k,1)
					endif
					if (k==20) then
						enorte=M(j,1,1)
					else
						enorte=M(j,k+1,1)
					endif
				else
					if (j==20) then
						eeste=M(1,k,l+1)
					else
						eeste=M(j+1,k,l+1)
					endif
					if (k==20) then
						enorte=M(j,1,l+1)
					else
						enorte=M(j,k+1,l+1)
					endif
				endif
				suma=norte+este+enorte+eeste
				E=E+2*suma
			else
				if (k==20) then
					norte=M(j,1,l)
				else
					norte=M(j,k+1,l)
				endif
				if (j==20) then
					este=M(1,k,l)
				else
					este=M(j+1,k,l)
				endif
				if (l==20) then
					if (j==20) then
						eeste=M(1,k,1)
					else
						eeste=M(j+1,k,1)
					endif
					if (k==20) then
						enorte=M(j,1,1)
					else
						enorte=M(j,k+1,1)
					endif
				else
					if (j==20) then
						eeste=M(1,k,l+1)
					else
						eeste=M(j+1,k,l+1)
					endif
					if (k==20) then
						enorte=M(j,1,l+1)
					else
						enorte=M(j,k+1,l+1)
					endif
				endif
				suma=norte+este+enorte+eeste
				E=E-2*suma
			endif
		enddo
	enddo
enddo

end subroutine Energia

subroutine datos(Etotal,E2total,E)

real (kind=dp), intent (inout) :: E,Etotal,E2total 

Etotal=Etotal+E
E2total=E2total+E*E

end subroutine datos

end module magnitudes
