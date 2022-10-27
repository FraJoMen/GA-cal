module Read_input
!DECLARATION OF VARIABLES
	implicit none

	integer, private :: i, j ! Counter for the do loops

	! Domain boundary of the search space
	real*8::phi_min,phi_max! min. and max. value of model parameter........φ
	real*8::hs_min,hs_max! min. and max. value of model parameter..........hs
	real*8::n_min,n_max! min. and max. value of model parameter............n
	real*8::ec0_min,ec0_max! min. and max. value of model parameter........ec0
	real*8::alpha_min,alpha_max! min. and max. value of model parameter....α
	real*8::beta_min,beta_max! min. and max. value of model parameter......β
	real*8::rat_ei_min,rat_ei_max! min. and max. value of model parameter..λi
	real*8::rat_ed_min,rat_ed_max! min. and max. value of model parameter..λd

	! GA set up
	real*8:: mu_in,mu_en ! Initial and final fraction of mutated individuals
	real*8:: N_F ! Number of individuals that are allowed to mate
	integer:: N_POP ! Population size
	integer:: P_dim ! Dimension of the optimization problem
	integer:: N_ITER,ITER ! Number of iteration and iteration counter
	integer:: N_E ! Number of elite individual passed to the next generation

	! Number of experimental data set and integration step
	integer:: Number_EDOTest,Number_TxDTest! number of eodometrics and triaxials
    !                                        drained test
	integer:: EDO_len,TxD_len! length of oedometric and triaxial input file
	integer:: init_EDO_len,init_TxD_len ! length of intial condition for the oe-
	!                                     dometric and triaxial input file
	integer:: N_Step ! numbers of steps for the integration

	! Initial condition
	real*8,allocatable:: init_EDO(:,:)! array of initial condition for eodomet-
	!                                   rics tests.This vector store:
	!                                   (T1 [kPa],T2 [kPa],e [-]).
	real*8,allocatable:: init_TxD(:,:) ! array of initial condition for triaxia-
	!                                   ls tests.This vector store:
	!                                   (T1 [kPa],T2 [kPa],e [-]).

	! Stored array variable
	real*8,allocatable:: X_EDO(:,:)! Array of array of experimental data
	!                                for the eodometrics tests. This vector sto-
	!                                re the evolution of state variable:
	!                                (T1 [MPa],T2 [MPa],e [-]).
    real*8,allocatable:: X_TxD(:,:)! Array of array of experimental data
	!                                for the triaxials tests. This vector store
	!                                the evolution of state variable:
	!                                (T1 [MPa],T2 [MPa],e [-]).
	! I/O variable
	character(len = 80):: File_title! Variable for skip the header of the file
	character(len = 30):: des,parameter_unit ! Variable for skip the description
	!                                          and unit of measure in the input
	!                                          file.
	character(len = 80):: data_Input! Name of file to read - Input file
	character(len = 80):: data_Domainfile! Name of file to read - Search space
	character(len = 80):: data_GAsetup! Name of file to read - GA set up
	character(len = 80):: data_EDO_init! Name of file to read - intials conditi-
	!                                    on for eodometrics tests.
	character(len = 80):: data_TxD_init! Name of file to read - intials conditi-
	!                                    on for triaxials tests.
	character(len = 80):: data_EDO,data_TxD! Name of file to read - data of
	!                                        eodometrics and triaxials tests.
	integer :: rdopst! Variable to store the iostat in the reading of file
	character(132):: line! Temporary variable to store a line from an input
	!                      file. Is used to get the length of file
	!                      data_EDO and data_TxD.

	! Evaluate Pop
	real*8::T1,T2! Principal and second Cauchy effective stress
	real*8::e! Void ratio
	real*8::f_d,f_s! Pyknotropy and barotropy coefficient
	real*8::e_d,e_c,e_i! The "Bauer's"  void ratio eq. (12)
	real*8::a! Coefficien compute by co_a function eq. (9)
	real*8::D1=-1! Principal component of the rate deformation tensors
    real*8::D2! Second component of the rate deformation tensors
	real*8::phi,hs,n,ec0,alpha,beta,ei0,ed0! Model parameters
	real*8::L11,L12,L21,L22! Linear part of the hypoplastic constitutive
	!                              matrix in axial symmetric condition
	real*8::N1,N2! Non linear part of the hypoplastic constitutive matrix in
	!              axial symmetric condition
	real*8::dydt_EDO(3)! gradient of the dinamic sistem for a oedometric test
    real*8::dydt_TxD(3)! gradient of the dinamic sistem for a triaxial drainded
    !                    test
	real*8::X_I,X_II! Possible solution of the quadratic equation ......
	real*8::tr_D,tr_T! Trace of Cauchy effective stress and rate deformation
	!                  tensors
    real*8::tr2_D,tr2_T! Trace of square of Cauchy effective stress and rate
	!                    deformation tensors
	real*8,allocatable::err_TxD_q(:)! Vector of errors on the Triaxial dev.
	!                                 plane(N_POP)
    real*8,allocatable::err_TxD_e(:)! Vector of errors on the Triaxial vol.
    !                                 plane(N_POP)
    real*8,allocatable::err_EDO(:)! Vector of errors on the oedometric plane
    !                               (N_POP)
	real*8,allocatable::Acurve_EDO(:,:,:,:)! array of the adimensional curve
	!                                    of the model prediction for triaxial
	!                                    deviatoric plane.
 	!                                    (N_Step,2,N_POP,Number_TxDTest).
	real*8,allocatable::Acurve_TxD_q(:,:,:,:)! array of the adimensional curve
	!                                    of the model prediction for triaxial
	!                                    deviatoric plane.
 	!                                    (N_Step,2,N_POP,Number_TxDTest).
    real*8,allocatable::Acurve_TxD_e(:,:,:,:)! array of the adimensional curve
	!                                    of the model prediction for triaxial
	!                                    volumetric plane.
 	!                                    (N_Step,2,N_POP,Number_TxDTest).
	real*8,allocatable::Apoint_TxD_q(:,:)
    real*8,allocatable::Apoint_TxD_e(:,:)
    real*8,allocatable::Apoint_EDO(:,:)! array of the adimensional curve
	!                                    of the model prediction for triaxial
	!                                    deviatoric plane.
 	!                                    (N_Step,2,N_POP,Number_TxDTest).
	real*8,allocatable:: mask_Point(:,:)!
	real*8,allocatable:: Xadm(:),Yadm(:)!
    real*8,allocatable::HX_EDO(:,:,:,:)! array (N_step,3,N_POP,Number_TxDTest)
    !                                    of the curve repose of the model predi-
    !                                    ction for eodometric tests. The numbers
    !                                    of colums 3 store the state vector:
	!                                    (T1 [MPa],T2 [MPa],e [-]).
    real*8,allocatable::HX_TxD(:,:,:,:)! array (N_step,3,N_POP,Number_TxDTest)
    !                                    of the curve repose of the model predi-
    !                                    ction for triaxial drained tests. The
    !                                    numbers of colums 3 store the state
    !                                    vector: (T1 [MPa],T2 [MPa],e [-]).
	real*8,allocatable::RX_TxD(:,:,:,:)! array (N_step,3,N_POP,Number_TxDTest)
	!                                    of the curve repose of the model predi-
	!                                    ction for triaxial drained test in
	!                                    terms of Roscoe variable.
	!                                    The numbers of colums 3 store the vec.:
	!                                    (q [MPa],epsilon_v [-],epsilon_a [-]).

	real*8,allocatable::Xadm_q(:),Yadm_q(:)! abscissa and ordinate for the
	!                                        dimensionlessness of the curves
	!                                        in the deviatoric plane.
	real*8,allocatable::Xadm_e(:),Yadm_e(:)! abscissa and ordinate for the
	!                                        dimensionlessness of the curves
	!                                        in the volumetric plane.
    real*8::Weight_EDO! Weight of errors in the oedometric plane
	real*8::Weight_TxD_q! Weight of errors in the Triaxial dev.
	real*8::Weight_TxD_e! Weight of errors in the Triaxial vol.
	real*8, allocatable::Score(:)! Vector for store the score of each individual
	!                              (N_POP)
	real*8,allocatable::time_TxD(:)! Vector whit the time integration for
	!                                triaxials tests (Number_TxDTest)
    real*8,allocatable::time_EDO(:)! Vector whit the time integration for
    !                                eodometrics tests (Number_EDOTest)

	integer::N_mask_Point
	integer:: row(1)
    integer:: mloc(1)
	integer,allocatable::Clas(:)

	logical,allocatable:: mk(:)
	logical, dimension(:),allocatable:: mask_EDO,mask_TxD
	logical,dimension(:,:),allocatable:: mask_TxD_dat
	logical, dimension(:,:),allocatable::mask_EDO_dat

	! Create Pop
	real*8, allocatable::POP(:,:)
	real*8, allocatable::POP_New(:,:)
	real*8, allocatable::Theta(:,:)
	real*8, allocatable::ranU1(:,:), ranT1(:,:)
	real*8, allocatable::ranU2(:,:), ranT2(:,:)
	real*8, allocatable::P_mat(:,:),P_ran(:,:),P_com(:,:),P_eli(:,:)
	real*8::L=1.0,M=1.0,mut
	real*8::Fc,U

	integer*8,dimension(:,:), allocatable::COMB_MAT
	integer*8,dimension(:), allocatable::X,Y
	integer:: N_mat,N_com,N_ran,N_mut


	contains
		subroutine Read_in()
! DISPLAY INITIAL HEADER
	 write(*,*)
	 write(*,*) "     _____                                      _ "
	 write(*,*) "    / ____|     /\                             | |"
     write(*,*) "   | |  __     /  \     ______    ___    __ _  | |"
	 write(*,*) "   | | |_ |   / /\ \   |______|  / __|  / _` | | |"
	 write(*,*) "   | |__| |  / ____ \           | (__  | (_| | | |"
	 write(*,*) "    \_____| /_/    \_\           \___|  \__,_| |_|"
	 write(*,*)
     write(*,*)


! READ THE INPUT FILE
	! prompt for input file name
	do
		write (*,'(a)', advance="no") "Input File Name string variable &
		                              & '...' >> 'path/filname.txt': "
		! read input file name
		read (*,*) data_Input
		! open input file (read access)
		! if open unsuccessful, display error message
		!  otherwise, end loop
	 open(10, file=data_Input, status="old", action="read", &
	      &position="rewind", iostat=rdopst )
	 if (rdopst==0) exit
	 write (*,'(a/,a)') "Unable to open input file.", "Please re-enter"
	 PAUSE
	 end do
	 write(*,*)
	!Read input information
	!Skip title
	read(10,*) File_title
	write(*,*) File_title
	write(*,*)
	!Extract value
	!Number_EDOTest
	read(10,*) des,Number_EDOTest
	write(*,*) des,Number_EDOTest
	!Number_TxDTest
	read(10,*) des,Number_TxDTest
	write(*,*) des,Number_TxDTest
	!Number of steps
	read(10,*) des,N_Step
	write(*,*) des,N_Step
	!Oedometric Wheigth
	read(10,*) des,Weight_EDO
	write(*,*) des,Weight_EDO
	!Triaxial deviatoric Wheigth
	read(10,*) des,Weight_TxD_q
	write(*,*) des,Weight_TxD_q
	!Triaxial deformation Wheigth
	read(10,*) des,Weight_TxD_e
	write(*,*) des,Weight_TxD_e
	!data_TxD,data/TxD_data.txt
	read(10,*) des,data_TxD
	write(*,*) des,data_TxD
	!data_EDO,data/EDO_data.txt
	read(10,*) des,data_EDO
	write(*,*) des,data_EDO
	!data_Domainfile,data/Domain.txt
	read(10,*) des,data_Domainfile
	write(*,*) des,data_Domainfile
	!data_GAsetup,data/GAsetup.txt
	read(10,*) des,data_GAsetup
	write(*,*) des,data_GAsetup
	!data_EDO_init,init_EDO_data.txt
	read(10,*) des,data_EDO_init
	write(*,*) des,data_EDO_init
	!data_TxD_init,init_TxD_data.txt
	read(10,*) des,data_TxD_init
	write(*,*) des,data_TxD_init
	write(*,*) '-------------------------------------------------------'
	write(*,*)
	close(10)

!READ PARAMETER'S DOMAIN FOR SAND HYPOPLASTICITY
	data_Domainfile='data/'// data_Domainfile
	  do
		! open input file (read access)
		! if open unsuccessful, display error message
		!  otherwise, end loop
		open(11, file=data_Domainfile, status="old", action="read",&
		     & position="rewind", iostat=rdopst )
		if (rdopst==0) exit
		write (*,'(a/,a)') "ERROR:: Unable to open data domain file.",&
		      & "Please check the input file."
        PAUSE
		STOP
	 end do
	! Skip title
	read(11,*) File_title
	write(*,*) File_title
	write(*,*)
	! Extract value
	!phi
	read(11,*) des,phi_min,phi_max,parameter_unit
	write(*,11) des,phi_min,phi_max,parameter_unit
	!hs
	read(11,*) des,hs_min,hs_max,parameter_unit
	write(*,11) des,hs_min,hs_max,parameter_unit
	!n
	read(11,*) des,n_min,n_max,parameter_unit
	write(*,11) des,n_min,n_max,parameter_unit
	!ec0
	read(11,*) des,ec0_min,ec0_max,parameter_unit
	write(*,11) des,ec0_min,ec0_max,parameter_unit
	!alpha
	read(11,*) des,alpha_min,alpha_max,parameter_unit
	write(*,11) des,alpha_min,alpha_max,parameter_unit
	!beta
	read(11,*) des,beta_min,beta_max,parameter_unit
	write(*,11) des,beta_min,beta_max,parameter_unit
	!rat_ed
	read(11,*) des,rat_ed_min,rat_ed_max,parameter_unit
	write(*,11) des,rat_ed_min,rat_ed_max,parameter_unit
	!rat_ei
	read(11,*) des,rat_ei_min,rat_ei_max,parameter_unit
	write(*,11) des,rat_ei_min,rat_ei_max,parameter_unit
	write(*,*) '-------------------------------------------------------'
	write(*,*)
	close(11)

!READ EODOMETRIC SPERIMENTAL DATA
	data_EDO='data/'// data_EDO
	do
		! if open unsuccessful, display error message
		!  otherwise, end loop
		open(12, file=data_EDO, status="old", action="read", &
		    & position="rewind", iostat=rdopst )
		if (rdopst==0) exit
		write (*,'(a/,a)') "ERROR:: Unable to open data domain file.",&
		      & "Please check the input file."
        PAUSE
		STOP
	 end do
	!Skip title
	read(12,*) File_title
	write(*,*) File_title
	read(12,*) File_title
	write(*,*) File_title
	write(*,*)
	!Extract value
	i = 1
	do
		! read line from input file
		read (12, '(a)', iostat=rdopst) line
		! if end of file, exit loop
		if (rdopst >0) stop "read error"
		if (rdopst <0) exit
		! count the length of data (length file -2) file
		i = i + 1
	 end do
	! X_EDO
	EDO_len=i-1
	allocate(X_EDO(EDO_len,3))
	rewind(12)
	read(12,*) File_title
	read(12,*) File_title
	i = 1
	do
		! read line from input file and alocate in a internal variable X_EDO
		read(12,*,iostat=rdopst) X_EDO(i,1),X_EDO(i,2),X_EDO(i,3)
		! if end of file, exit loop
		if (rdopst >0) stop "read error"
		if (rdopst <0) exit
		i = i + 1
	 end do
	!Print the array X_EDO
	do i=1,EDO_len
		write(*,*) i, X_EDO(i,:)
	enddo
	write(*,*) '-------------------------------------------------------'
	write(*,*)
	close(12)

!READ TRIAXIAL DRAINED SPERIMENTAL DATA
	data_TxD='data/'// data_TxD
	do
		! open input file (read access)
		! if open unsuccessful, display error message
		! otherwise, end loop
		open(13, file=data_TxD, status="old", action="read", &
		     &position="rewind", iostat=rdopst )
		if (rdopst==0) exit
		write (*,'(a/,a)') "ERROR:: Unable to open data domain file.",&
		                   & "Please check the input file."
        PAUSE
		STOP
	end do
	!Skip title
	read(13,*) File_title
	write(*,*) File_title
	read(13,*) File_title
	write(*,*) File_title
	write(*,*)
	!Extract value
	i = 1
	do
		! read line from input file
		read (13, '(a)', iostat=rdopst) line
		! if end of file, exit loop
		if (rdopst >0) stop "read error"
		if (rdopst <0) exit
		! count the length of data (length file -2) file
		i = i + 1
	end do
	! X_TxD
	TxD_len=i-1
	allocate(X_TxD(TxD_len,4))
	rewind(13)
	read(13,*) File_title
	read(13,*) File_title
	i = 1
	do
		! read line from input file and alocate in a internal variable X_TxD
		read(13,*,iostat=rdopst) X_TxD(i,1),X_TxD(i,2),X_TxD(i,3),X_TxD(i,4)
		! if end of file, exit loop
		if (rdopst >0) stop "read error"
		if (rdopst <0) exit
		i = i + 1
	end do
	!Print the array X_EDO
	do i=1,TxD_len
		write(*,*) i, X_TxD(i,:)
	end do
	write(*,*) '-------------------------------------------------'
	write(*,*)
	close(13)

!READ INITIAL CONDITION EDO TESTS
	data_EDO_init='data/'// data_EDO_init
	do
		! open input file (read access)
		! if open unsuccessful, display error message
		!  otherwise, end loop
		open(14, file=data_EDO_init, status="old", action="read", &
		     &position="rewind", iostat=rdopst )
		if (rdopst==0) exit
		write (*,'(a/,a)') "ERROR:: Unable to open data domain file.",&
		                   & "Please check the input file."
        PAUSE
		STOP
	 end do
	!Skip title
	read(14,*) File_title
	write(*,*) File_title
	read(14,*) File_title
	write(*,*) File_title
	write(*,*)
	!Extract value
	i = 1
	do
		! read line from input file
		read (14, '(a)', iostat=rdopst) line
		! if end of file, exit loop
		if (rdopst >0) stop "read error"
		if (rdopst <0) exit
		! count the length of data
		i = i + 1
	end do
	! init_EDO
	init_EDO_len=i-1
	allocate(init_EDO(init_EDO_len,3))
	rewind(14)
	read(14,*) File_title
	read(14,*) File_title
	i = 1
	do
		! read line from input file and alocate
		! in a internal variable X_TxD
		read(14,*,iostat=rdopst) init_EDO(i,1),init_EDO(i,2),init_EDO(i,3)
		! if end of file, exit loop
		if (rdopst >0) stop "read error"
		if (rdopst <0) exit
		i = i + 1
	end do
	!Print the array init_EDO
	do i=1,init_EDO_len
		write(*,*) i, init_EDO(i,:)
	end do
	write(*,*) '-------------------------------------------------'
	write(*,*)
	close(14)

!READ INITIAL CONDITION TXD TESTS
	data_TxD_init='data/'// data_TxD_init
	do
		! open input file (read access)
		! if open unsuccessful, display error message
		!  otherwise, end loop
		open(15, file=data_TxD_init, status="old", action="read", &
		     &position="rewind", iostat=rdopst )
		if (rdopst==0) exit
		write (*,'(a/,a)') "ERROR:: Unable to open data domain file.",&
		                   & "Please check the input file."
        PAUSE
		STOP
	end do
	!Skip title
	read(15,*) File_title
	write(*,*) File_title
	read(15,*) File_title
	write(*,*) File_title
	write(*,*)
	!Extract value
	i = 1
	do
		! read line from input file
		read (15, '(a)', iostat=rdopst) line
		! if end of file, exit loop
		if (rdopst >0) stop "read error"
		if (rdopst <0) exit
		! count the length of data
		i = i + 1
	end do
	! init_TxD
	init_TxD_len=i-1
	allocate(init_TxD(init_TxD_len,3))
	rewind(15)
	read(15,*) File_title
	read(15,*) File_title
	i = 1
	do
		! read line from input file and alocate
		! in a internal variable X_TxD
		read(15,*,iostat=rdopst) init_TxD(i,1),init_TxD(i,2),init_TxD(i,3)
		! if end of file, exit loop
		if (rdopst >0) stop "read error"
		if (rdopst <0) exit
		i = i + 1
	end do
	! Print the array init_EDO
	do i=1,init_TxD_len
		write(*,*) i, init_TxD(i,:)
	end do
	write(*,*) '-------------------------------------------------'
	write(*,*)
	close(15)

!READ GA SET UP
	data_GAsetup='data/'// data_GAsetup
	do
		! open input file (read access)
		! if open unsuccessful, display error message
		! otherwise, end loop
		open(14, file=data_GAsetup, status="old", action="read",&
		    & position="rewind", iostat=rdopst )
		if (rdopst==0) exit
		write (*,'(a/,a)') "ERROR:: Unable to open GA set up.", &
		                 & "Please check the input file."
        PAUSE
		STOP
	end do
	! Skip title
	read(14,*) File_title
	write(*,*) File_title
	write(*,*)
	! Extract value
	read(14,*) des,N_POP
	write(*,14) des,N_POP
	read(14,*) des,P_dim
	write(*,14) des,P_dim
	read(14,*) des,N_ITER
	write(*,14) des,N_ITER
	read(14,*) des,N_E
	write(*,14) des,N_E
	read(14,*) des,N_F
	write(*,15) des,N_F
	read(14,*) des,mu_en
	write(*,15) des,mu_en
	read(14,*) des,mu_in
	write(*,15) des,mu_in
	write(*,*) '------------------------------------------'
	write(*,*)
	close(14)

! FORMAT
	 11 FORMAT ( A,f10.2,f10.2,3x,A )
	 14 FORMAT ( A,5x,i5.1)
	 15 FORMAT ( A,3x,f10.2)

! COMPUTE THE INTEGRATION TIME
	!Oedometric test
	allocate(time_EDO(Number_EDOTest))
	allocate(mask_EDO(size(X_EDO)/3))
    do i =1,Number_EDOTest
		mask_EDO=.FALSE.
		where (X_EDO(:,3)==real(i))
			mask_EDO=.TRUE.
		end where
	time_EDO(i)=-log((1+MINVAL(X_EDO(:,1),mask_EDO))/(1+init_EDO(i,3)))&
	           &*1.01
	end do
 	!Triaxial test
	allocate(time_TxD(Number_TxDTest))
	allocate(mask_TxD(size(X_TxD)/4))
    do i =1,Number_TxDTest
		mask_TxD=.FALSE.
		where (X_TxD(:,4)==real(i))
			mask_TxD=.TRUE.
		end where
	time_TxD(i)=MAXVAL(X_TxD(:,1),mask_TxD)/100.0*1.05
    end do

!! ALOCATION OF VARIABLE
	! Evaluate Pop
	allocate(Score(N_POP))
	allocate(Clas(N_POP))
	allocate(err_EDO(N_POP))
	allocate(HX_EDO(N_step,3,N_POP,Number_EDOTest))
	allocate(HX_TxD(N_step,3,N_POP,Number_TxDTest))
	allocate(RX_TxD(N_step,3,N_POP,Number_TxDTest))
	allocate(Acurve_EDO(N_Step,2,N_POP,Number_EDOTest))
	allocate(Apoint_EDO(EDO_len,2))
	allocate(mask_EDO_dat(EDO_len,Number_EDOTest))
    allocate(Xadm(Number_EDOTest))
    allocate(Yadm(Number_EDOTest))
    allocate(err_TxD_q(N_POP))
    allocate(err_TxD_e(N_POP))
    allocate(Acurve_TxD_q(N_Step,2,N_POP,Number_TxDTest))
    allocate(Acurve_TxD_e(N_Step,2,N_POP,Number_TxDTest))
    allocate(Apoint_TxD_q(TxD_len,2))
    allocate(Apoint_TxD_e(TxD_len,2))
    allocate(mask_TxD_dat(TxD_len,Number_TxDTest))
    allocate(Xadm_q(Number_TxDTest))
    allocate(Yadm_q(Number_TxDTest))
    allocate(Xadm_e(Number_TxDTest))
    allocate(Yadm_e(Number_TxDTest))
    allocate(mk(N_POP))
	! Create Pop
	allocate(POP(N_POP,P_dim))
	allocate(POP_New(N_POP,P_dim))

		end subroutine Read_in
end module Read_input
