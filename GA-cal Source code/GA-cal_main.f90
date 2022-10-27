program main
	use Read_input
	use Create_Pop
	use Evaluate_Pop
	use Write_out
    real :: start, finish
    call cpu_time(start)
	100 FORMAT ( 5x,A,8x,A,8x,A,9x,A,7x,A,5x,A,6x,A,4x,A)
	101 FORMAT ( f10.2,f10.2,f10.2,f10.2,f10.2,f10.2,f10.2,f10.2)
	200 FORMAT ( 4x,A,3x,A,8x,A)
	201 FORMAT ( i8,ES14.7,i8)
! READ INPUT DATA
	call Read_in
! CREATE AND EVALUATE A INITIAL POPULATION
    call Initial_Pop
	print*, 'Evaluate the initial Pop:'
	call intE_EDO
	call intE_TxD
	call elab_TxD
	call COST_EDO
	call COST_TxD
	call Cost_Pop
	call merge_argsort
	write(*,*)
! GA OPTIMIZATION
	do ITER=1,N_ITER
		write(*,*)
		print*, '------------------------------------------------------'
		print*, '- ITERATION ----------------',ITER,'/',N_ITER
		write(*,*)
		call Update_Pop
		! Evaluate and update population
		POP=POP_New
		call intE_EDO
		call intE_TxD
		call elab_TxD
		call COST_EDO
		call COST_TxD
		call Cost_Pop
		call merge_argsort
	end do
	call merge_argsort
! PRINT THE SCORE AND PARAMETERS VALUE - FINAL POPULATION
	write(*,*)
	write(*,*) '#-----------Evaluated Pop-------------------------------'
	write(*,200)  'Class', 'Score  ', 'ID'
	write(*,*) '-------------------------------------------------------'
	do i=1,10
	    write(*,201)  i,Score(Clas(i)),Clas(i)
	end do
	write(*,*) '----------------------------------'
	write(*,*)
	write(*,*) '#-----------Updated Pop -------------------------------'
	write(*,100) 'phi','hs','n','ec0','alpha','beta','rat_ei','rat_ed'
	write(*,*) '-------------------------------------------------------'
	do i=1,6
		write(*,101) POP(Clas(i),:)
	end do
	write(*,*)
	write(*,*) '#-----------Best fit-----------------------------------'
	write(*,100) 'phi','hs','n','ec0','alpha','beta','rat_ei','rat_ed'
	write(*,*) '-------------------------------------------------------'
	write(*,101) POP(Clas(1),:)
! WTRITE THE OUTPUT
	call Write_ou
	call cpu_time(finish)
	write(*,*)
    write(*,*) 'CPU Time =',finish-start,' seconds.'
    write(*,*) '... Press Enter key to quit the program ...'
    read(*,*)
end program main
