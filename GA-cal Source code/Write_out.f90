
module Write_out
!Declaration of variables
	use Read_input
	contains

! WRITE THE OUTPUT
		subroutine Write_ou
	!write(*,*) X_EDO(1,1) sream = all information are sequential
	open(50, file="data_structure.dat", access="stream")
	write(50) EDO_len, TxD_len,N_Step,N_POP,Number_EDOTest,Number_TxDTest
	close(50)
	!write(*,*) X_EDO(1,1) sream = all information are sequential
	open(51, file="X_EDO.dat", access="stream")
	write(51) X_EDO
	close(51)
	!write(*,*) X_TxD(1,1) sream = all information are sequential
	open(52, file="X_TxD.dat", access="stream")
	write(52) X_TxD
	close(52)
	!write(*,*) HX_EDO sream = all information are sequential
	open(53, file="HX_EDO.dat", access="stream")
	write(53) HX_EDO
	close(53)
	!write(*,*) HX_TxD sream = all information are sequential
	open(53, file="RX_TxD.dat", access="stream")
	write(53) RX_TxD
	close(53)

    100 FORMAT ( 5x,A,8x,A,8x,A,9x,A,7x,A,5x,A,6x,A,4x,A)
	101 FORMAT ( f10.2,f10.2,f10.2,f10.2,f10.2,f10.2,f10.2,f10.2)
    open(unit=100,file='best_fit.txt',form='formatted',&
	   & status='unknown')
	write(100,*) '#-----------Best fit-----------------------------------'
	write(100,100) 'phi','hs','n','ec0','alpha','beta','rat_ei','rat_ed'
	write(*,*) '-------------------------------------------------------'
	write(100,101) POP(Clas(1),:)
    close(100)

		end subroutine Write_ou
end module Write_out
