module Create_Pop
! DECLARATION OF VARIABLES
    use Read_input
	integer,private:: values(1:8),kkk ! variables to set a different seed  via clock
	integer,private::i,j,k
	integer,dimension(:),allocatable,private:: seed
	real*8,dimension(:,:),allocatable,private::Omega ! Matrix of random numbers used to compute the initial population and the re-initialized population.
	real*8,private::Lower,Upper

	contains
! DEFINE THE INITIAL POPULATION
		subroutine Initial_Pop
    !DESCRIPTION:
	!Compute the initial population by continuous uniform distribution between
	!the Lower and Upper values of each model parameter defined in the input.

	100 FORMAT ( 5x,A,8x,A,8x,A,9x,A,7x,A,5x,A,6x,A,4x,A)
	101 FORMAT ( f10.2,f10.2,f10.2,f10.2,f10.2,f10.2,f10.2,f10.2)

	! Compute a random vector
	call date_and_time(values=values)
	call random_seed(size=kkk) !set the seed
	allocate(seed(1:kkk))
	seed(:)=values(8)*values(7)+values(6)*values(5)
	! Random uniform
	call random_seed(put=seed)
	allocate(Omega(N_POP,P_dim))
	call random_number(Omega)
	! 1 -phi
	Upper=phi_max
	Lower=phi_min
	do i=1,N_POP
		POP(i,1)=(Upper-Lower)*Omega(i,1)+Lower
	end do
	! 2 -hs
	Upper=hs_max
	Lower=hs_min
	do i=1,N_POP
		POP(i,2)=(Upper-Lower)*Omega(i,2)+Lower
	end do
	! 3 -n
	Upper=n_max
	Lower=n_min
	do i=1,N_POP
		POP(i,3)=(Upper-Lower)*Omega(i,3)+Lower
	end do
	! 4 -ec0
	Upper=ec0_max
	Lower=ec0_min
	do i=1,N_POP
		POP(i,4)=(Upper-Lower)*Omega(i,4)+Lower
	end do
	! 5 -alpha
	Upper=alpha_max
	Lower=alpha_min
	do i=1,N_POP
		POP(i,5)=(Upper-Lower)*Omega(i,5)+Lower
	end do
	! 6 -beta
	Upper=beta_max
	Lower=beta_min
	do i=1,N_POP
		POP(i,6)=(Upper-Lower)*Omega(i,6)+Lower
	end do
	! 7 -rat_ei
	Upper=rat_ei_max
	Lower=rat_ei_min
	do i=1,N_POP
		POP(i,7)=(Upper-Lower)*Omega(i,7)+Lower
	end do
	! 8 -rat_ed
	Upper=rat_ed_max
	Lower=rat_ed_min
	do i=1,N_POP
		POP(i,8)=(Upper-Lower)*Omega(i,8)+Lower
	end do
	deallocate(seed)
	deallocate(Omega)
	open(unit=100,file='log_PoP.txt',form='formatted',&
	   & status='unknown')
	write(100,*) '#--- Initial_population -----------------------------'
	write(100,100) 'Omega','hs','n','ec0','alpha','beta','rat_ed','rat_ei'
	write(100,*) '-----------------------------------------------------'
	do i=1,N_POP
	    write(100,101) POP(i,:)
	end do
	write(100,*) '-----------------------------------------------------'
	close(100)
			end subroutine Initial_Pop
! DEFINE THE NEW POPULATION
		subroutine Update_Pop
    !DESCRIPTION:
	!Combines elitism, mutation, selection and cross-over to generate and define
	!a new population with an improved average cost.

	100 FORMAT ( 5x,A,8x,A,8x,A,9x,A,7x,A,5x,A,6x,A,4x,A)
	101 FORMAT ( f10.2,f10.2,f10.2,f10.2,f10.2,f10.2,f10.2,f10.2)
	104 FORMAT ( f10.2,f10.2)
	200 FORMAT ( 9x,A,9x,A,8x,A,2x,A)
	201 FORMAT ( i10,i10,i10,i10)
	301 FORMAT ( i10,i10,i10,i10,i10,i10,i10,i10)

	open(unit=100,file='log_PoP.txt',form='formatted',&
	   & status='unknown',access='append')

	! Compute the Mutation ratio
	U=N_F*N_POP
	mut= mu_in*(mu_en/mu_in)**((real(ITER)-1.0)/(real(N_ITER)-1.0))
	N_mut=int((mut)*real(N_POP)) !
	N_com=N_mut*real(ITER)/real(N_ITER)
	N_ran=N_mut-N_com
	N_mat=N_POP-N_E-N_com-N_ran

	! Alocate the vaiable
    allocate(P_mat(N_mat,P_dim))
    allocate(P_ran(N_ran,P_dim))
    allocate(P_com(N_com,P_dim))
    allocate(P_eli(N_E,P_dim))

    ! Compute elite Population
    do i=1,N_E
		do k=1,P_dim
        P_eli(i,k)=POP(Clas(i),k)
	    end do
	end do
	write(100,*) '#---------Update Pop P_eli-----',ITER,'/',N_ITER
	write(100,100) 'Omega','hs','n','ec0','alpha','beta','rat_ei','rat_ed'
	write(100,*) '-----------------------------------------------------'
	do i=1,N_E
	    write(100,101) P_eli(i,:)
	end do
	write(100,*)
	! Compute mate Population
	call date_and_time(values=values)
	call random_seed(size=kkk)
	allocate(seed(1:kkk))
	seed(:)=values(8)*values(7)+values(6)*values(5)
	allocate(Theta(N_mat,P_dim))
	call random_number(Theta)
	! Vector of selected individual
	allocate(X(N_mat))
	allocate(Y(N_mat))
	allocate(ranU1(N_mat,2))
	call random_number(ranU1)
	! Triangular random distribution
	allocate(ranT1(N_mat,2))
	Fc=(M-L)/(U-L)
	do i=1,N_mat
	  do j=1,2
	   if( ranU1(i,j) < Fc) then
	   ranT1(i,j)= L+sqrt(ranU1(i,j)*(U-L)*(M-L))
	   else
	   ranT1(i,j)= U-sqrt((1-ranU1(i,j))*(U-L)*(U-M))
	   end if
	  end do
	end do
	X=int(ranT1(:,1))
	Y=int(ranT1(:,2))
	! Mask zero value
    where(X==0)
		X=1
    end where
    where(Y==0)
		Y=1
    end where
    ! Compute the population mated
	do i=1,N_mat
		do k=1,P_dim
        P_mat(i,k)=POP(Clas(X(i)),k)*Theta(i,k)+&
                   (1-Theta(i,k))*POP(Clas(Y(i)),k)
	    end do
	end do
	write(100,*)
	write(100,*) '#---------Update Pop P_mat-----',ITER,'/',N_ITER
	write(100,100) 'Omega','hs','n','ec0','alpha','beta','rat_ei','rat_ed'
	write(100,*) '-----------------------------------------------------'
	do i=1,N_mat
	    write(100,101) P_mat(i,:)
	end do
	write(100,*)

	! Create mutated Population
	seed(:)=values(8)*values(7)+values(6)*values(5)
	call random_seed(put=seed)
	allocate(Omega(N_ran,P_dim))
	call random_number(Omega)
	! 1 - phi
	Upper=phi_max
	Lower=phi_min
	do i=1,N_ran
		P_ran(i,1)=(Upper-Lower)*Omega(i,1)+Lower
	end do
	! 2 -hs
	Upper=hs_max
	Lower=hs_min
	do i=1,N_ran
		P_ran(i,2)=(Upper-Lower)*Omega(i,2)+Lower
	end do
	! 3 -n
	Upper=n_max
	Lower=n_min
	do i=1,N_ran
		P_ran(i,3)=(Upper-Lower)*Omega(i,3)+Lower
	end do
	! 4 -ec0
	Upper=ec0_max
	Lower=ec0_min
	do i=1,N_ran
		P_ran(i,4)=(Upper-Lower)*Omega(i,4)+Lower
	end do
	! 5 -alpha
	Upper=alpha_max
	Lower=alpha_min
	do i=1,N_ran
		P_ran(i,5)=(Upper-Lower)*Omega(i,5)+Lower
	end do
	! 6 -beta
	Upper=beta_max
	Lower=beta_min
	do i=1,N_ran
		P_ran(i,6)=(Upper-Lower)*Omega(i,6)+Lower
	end do
	! 7 -rat_ei
	Upper=rat_ei_max
	Lower=rat_ei_min
	do i=1,N_ran
		P_ran(i,7)=(Upper-Lower)*Omega(i,7)+Lower
	end do
	! 8 -rat_ed
	Upper=rat_ed_max
	Lower=rat_ed_min
	do i=1,N_ran
		P_ran(i,8)=(Upper-Lower)*Omega(i,8)+Lower
	end do
	write(100,*) '#-----------Update Pop P_ran-----',ITER,'/',N_ITER
	write(100,100) 'phi','hs','n','ec0','alpha','beta','rat_ei','rat_ed'
	write(100,*) '-------------------------------------------------------'
	do i=1,N_ran
	    write(100,101) P_ran(i,:)
	end do
	write(100,*)

	! Create mutated combination population
	allocate(COMB_MAT(N_com,P_dim))
	allocate(ranU2(N_com,P_dim))
	allocate(ranT2(N_com,P_dim))
	call random_number(ranU2)
	! Triangular random distribution
	U=N_POP
	Fc=(M-L)/(U-L)
	do i=1,N_com
	  do j=1,P_dim
	   if( ranU2(i,j) < Fc) then
	   ranT2(i,j)= L+sqrt(ranU2(i,j)*(U-L)*(M-L))
	   else
	   ranT2(i,j)= U-sqrt((1-ranU2(i,j))*(U-L)*(U-M))
	   end if
	  end do
	end do
	COMB_MAT=int(ranT2)
	! Mask zero value
    where(COMB_MAT==0)
		COMB_MAT=1
    end where
    ! Create P_com
    do i=1,N_com
		do k=1,P_dim
        P_com(i,k)=POP(Clas(COMB_MAT(i,k)),k)
	    end do
	end do
	write(100,*) '#-----------Update Pop P_com-----',ITER,'/',N_ITER
	write(100,100) 'Omega','hs','n','ec0','alpha','beta','rat_ed','rat_ei'
	write(100,*) '-------------------------------------------------------'
	do i=1,N_com
	    write(100,101) P_com(i,:)
	end do
	write(100,*)

	! Define updated population
	POP_New(1:N_E,:)=P_eli
	POP_New(1+N_E:1+N_E+N_mat,:)=P_mat
	POP_New(1+N_E+N_mat:1+N_E+N_mat+N_com,:)=P_com
	POP_New(1+N_E+N_mat+N_com:1+N_E+N_mat+N_com+N_ran,:)=P_ran
    ! Dealocate local variable
	deallocate(Theta)
	deallocate(COMB_MAT)
	deallocate(X)
	deallocate(Y)
	! Random numbers
	deallocate(seed)
	! Mate Pop
	deallocate(ranU1)
	deallocate(ranT1)
	deallocate(ranU2)
	deallocate(ranT2)
	deallocate(P_mat)
	deallocate(P_ran)
	deallocate(P_com)
	deallocate(P_eli)
	! Mutation Pop
	deallocate(Omega)

        end subroutine  Update_Pop
end module Create_Pop
