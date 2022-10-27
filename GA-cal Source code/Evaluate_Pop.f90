module Evaluate_Pop
! DECLARATION OF VARIABLES
	use Read_input
	integer,private:: i,j,k
	real*8, parameter :: pi=4.0*ATAN(1.0)
	contains

! COST FUNCTION
		subroutine Cost_Pop
    !DESCRIPTION:
	!It combines the computed errors with the SUBROUTINES COST_EDO and
	!COST_TxD for the definition of the score assigned to each
	!individual of the population. The higher it is, the greater the
	!deviation of the individual from the experimental data
	!
	Score=0
	do i=1,N_POP
	    Score(i)=Weight_EDO*err_EDO(i)    + &
	             Weight_TxD_q*err_TxD_q(i)+ &
	             Weight_TxD_e*err_TxD_e(i)
	end do
	print*, 'min - max delta OE  ', MINVAL(err_EDO),MAXVAL(err_EDO)
	print*, 'min - max delta TD_q', MINVAL(err_TxD_q),MAXVAL(err_TxD_q)
	print*, 'min - max delta TD_e', MINVAL(err_TxD_e),MAXVAL(err_TxD_e)
	print*,'----------------------------'
	print*,'best score: ',MINVAL(Score)
		end subroutine Cost_Pop

! ARG SHOTR OF A ARRAY VECTOR
		subroutine merge_argsort
    !DESCRIPTION:
	!Perform a sort of the Score of the population. Returns an array of indices
	!that sorts the Score of the population in ascending order.

	mk = .TRUE.
	do i = 1, N_POP
		mloc=MINLOC(Score,mk)
		Clas(i)=mloc(1)
		mk(mloc) = .FALSE.
	end do
		end subroutine  merge_argsort

! CALCULATION OF THE TRACE
		function tr(AA,BB)
	!DESCRIPTION:
	!Compute the trace of a matrix M
	!      |A 0 0|
	! M= : |0 B 0|
	!      |0 0 B|
	implicit none
	real*8::AA,BB
	real*8::tr
	tr = AA+2*BB
    	end function tr

! CALCULATION OF THE SQUARE TRACE
    	function tr2(AA,BB)
    !DESCRIPTION:
	!Compute the trace of the square of the matrix M (M^2)
	!      |A 0 0|
	! M= : |0 B 0|
	!      |0 0 B|
	implicit none
	real*8::AA,BB
	real*8::tr2
	tr2 = AA*AA+2*BB*BB
    	end function tr2

! CALCULATION OF THE COEFFICIENT a
    	function co_a(phi)
	!DESCRIPTION:
	!Compute the coefficient a of SH model from eq.......
	real*8::phi
	real*8::co_a
	co_a = sqrt(3.0)*(3-sin(phi*pi/180))/(2*sqrt(2.0)*sin(phi*pi/180))
    	end function co_a

! BAWER'S LAWS
		function e_dci(e_dci0,tr_T,hs,n)
	!DESCRIPTION:
	!Compute the Bauer’s law for determinate the coefficient e_i, e_d,
	!and e_c of the SH model. See equation ....
	real*8::e_dci0,tr_T,hs,n
	real*8::e_dci
	e_dci=e_dci0*exp(-(-tr_T/(hs*1000))**n)
    	end function e_dci

! PYKNOTROPY COEFFICIENT f d
		function fd(tr_T,e,ed,ec,alpha)
	!DESCRIPTION:
	!Compute the pyknotropy coefficient f_d in the equation.....
	real*8::tr_T,e,ed,ec,alpha
	real*8::fd
	fd=((e-ed)/(ec-ed))**alpha
		end function fd

! BAROTROPY COEFFICIENT f_s
		function fs(tr_T,e,ed,ec,ei,ed0,ec0,ei0,a,hs,n,alpha,betha)
	!DESCRIPTION:
	!Compute the barotropy coefficient f_s in the equation.....
	real*8::tr_T,e,ed,ec,ei,ed0,ec0,ei0,a,hs,n,alpha,betha
	real*8::fs
	fs=(hs*1000)/n*(ei/e)**betha*(1+ei)/ei*(-tr_T/(hs*1000))&
	    &**(1-n)/(3+a**2-a*sqrt(3.0)*((ei0-ed0)/(ec0-ed0))**alpha)
		end function fs

! FRECHET DISTANT
		function Frechet_dist(N_Step,N_mask_Point,Point,Curve)
	!DESCRIPTION:
	!Computed the root mean square of the Fréchet distances of a
	!<<Points>> and <<Curve>>.
	integer::N_Step,N_mask_Point
	real*8,dimension(N_mask_Point,2)::Point
	real*8,dimension(N_Step,2)::Curve
	real*8::Frechet_dist
	real*8::dist(N_Step)
	real*8::Frechet(N_mask_Point)
	real*8::q,m
	logical, dimension(N_Step):: mk
	integer,dimension(2):: ID_min
	integer,dimension(1):: row
	integer::ii,kk,jj

	do kk=1,N_mask_Point
		do ii=1,N_Step
			dist(ii)=sqrt((Curve(ii,1)-Point(kk,1))**2+&
			            &(Curve(ii,2)-Point(kk,2))**2)
		end do
	! Find nearest points
		mk = .TRUE.
		do jj = 1, 2
			row=MINLOC(dist,mk)
			ID_min(jj)=row(1)
			mk(row) = .FALSE.
		end do
	! Frechet computation
		if(Curve(ID_min(1),1)==Curve(ID_min(2),1)) then
			Frechet(k)=abs(Curve(ID_min(1),1)-Point(k,1))
		else
		m=(Curve(ID_min(2),2)-Curve(ID_min(1),2))/&
		  (Curve(ID_min(2),1)-Curve(ID_min(1),1))
		q=(Curve(ID_min(2),1)*Curve(ID_min(1),2)-&
		  Curve(ID_min(1),1)*Curve(ID_min(2),2))&
		  /(Curve(ID_min(2),1)-Curve(ID_min(1),1))
		Frechet(kk)=abs(Point(kk,2)-(m*Point(kk,1)+q))/sqrt(1+m**2)
		end if
	end do

	Frechet_dist=norm2(Frechet)/N_mask_Point
		end function Frechet_dist

! EODOMETRIC GRADIENTE
		subroutine dydtEDO
    !DESCRIPTION:
	!Compute the gradient of oedometric test. See eq. .....
	!
	tr_T=tr(T1,T2)
    tr2_T=tr2(T1,T2)
    L11=tr_T**2/tr2_T*(1+a**2*T1**2/tr_T**2)
    L21=(a**2*T1*T2)/tr2_T
    N1=tr_T**2/tr2_T*(a/3*(5*T1-2*T2)/tr_T)
    N2=tr_T**2/tr2_T*(a/3*(4*T2-T1)/tr_T)
	dydt_EDO(1)= f_s*L11*D1+f_s*f_d*N1*sqrt(D1**2)
    dydt_EDO(2)= f_s*L21*D1+f_s*f_d*N2*sqrt(D1**2)
    dydt_EDO(3)= (1+e)*D1
		end subroutine dydtEDO

! EULER-EXPLICIT INTEGRATION OF THE OEDOMETRIC TEST
		subroutine intE_EDO
	! Perform the explicit Euler integration to simulate oedometers test using
	! dydtEDO subroutine

	HX_EDO=0
	do k=1,Number_EDOTest
	Dtime=time_EDO(k)/N_Step
		do i=1,N_POP
		! Model parameters
		phi=POP(i,1)
		hs =POP(i,2)
		n  =POP(i,3)
		ec0=POP(i,4)
		alpha=POP(i,5)
		beta=POP(i,6)
		ei0=POP(i,7)*ec0
		ed0=POP(i,8)*ec0
		! State variable
		T1=-init_EDO(k,1)/1000
		T2=-init_EDO(k,2)/1000
		e = init_EDO(k,3)
		tr_T=tr(T1,T2)
		e_d=e_dci(ed0,tr_T,hs,n)
		e_c=e_dci(ec0,tr_T,hs,n)
		e_i=e_dci(ei0,tr_T,hs,n)
		a=co_a(phi)
		f_d=fd(tr_T,e,e_d,e_c,alpha)
		f_s=fs(tr_T,e,e_d,e_c,e_i,ed0,ec0,ei0,a,hs,n,alpha,beta)
			do j=1,N_Step
			HX_EDO(j,1,i,k)=T1
			HX_EDO(j,2,i,k)=T2
			HX_EDO(j,3,i,k)=e
			! Check - tr_T < 0 e ei<e<ed
			if (tr_T >0) exit
			if (e>e_i) exit
			if (e<e_d) exit
			call dydtEDO
			T1=dydt_EDO(1)*Dtime+T1
			T2=dydt_EDO(2)*Dtime+T2
			e =dydt_EDO(3)*Dtime+e
			tr_T=tr(T1,T2)
			e_d=e_dci(ed0,tr_T,hs,n)
			e_c=e_dci(ec0,tr_T,hs,n)
			e_i=e_dci(ei0,tr_T,hs,n)
			a=co_a(phi)
			f_d=fd(tr_T,e,e_d,e_c,alpha)
			f_s=fs(tr_T,e,e_d,e_c,e_i,ed0,ec0,ei0,a,hs,n,alpha,beta)
			end do
		end do
	end do
		end subroutine intE_EDO

! TRIAXIAL GRADIENT
		subroutine dydtTxD
	!DESCRIPTION:
	!Compute the gradient of triaxial drained test. See eq. .....
	!

	tr_T=tr(T1,T2)
    tr2_T=tr2(T1,T2)
    L11=tr_T**2/tr2_T*(1+a**2*T1**2/tr_T**2)
    L12=2*a**2*T1*T2/tr2_T
    L21=(a**2*T1*T2)/tr2_T
    L22=tr_T**2/tr2_T*(1+a**2*2*T2**2/tr_T**2)
    N1=tr_T**2/tr2_T*(a/3*(5*T1-2*T2)/tr_T)
    N2=tr_T**2/tr2_T*(a/3*(4*T2-T1)/tr_T)
    ! Solve the quadratic equation
    X_I=-(L22*sqrt(-2*N2**2*D1**2*f_d**2+L22**2*D1**2+2*D1**2*L21**2)+2*D1*N2*L21*f_d)/(2*N2**2*f_d**2-L22**2)
    X_II=(L22*sqrt(-2*N2**2*D1**2*f_d**2+L22**2*D1**2+2*D1**2*L21**2)-2*D1*N2*L21*f_d)/(2*N2**2*f_d**2-L22**2)
    if (X_I > 0 .AND. X_II < 0) then
		D2= -(N2*f_d*X_I+D1*L21)/(L22)
	else if (X_II>0 .AND. X_I < 0) then
		D2 = -(N2*f_d*X_II+D1*L21)/(L22)
	else
		print*,'Not unique solution'
		D2 = 999
	end if
    ! Compute the gradient
    if (D2 == 999) then
		dydt_TxD(1)= 0
		dydt_TxD(2)= 0
		dydt_TxD(3)= 0
	else
		tr_D =tr(D1,D2)
		tr2_D=tr2(D1,D2)
		dydt_TxD(1)= f_s*L11*D1+f_s*L12*D2+f_s*f_d*N1*sqrt(tr2_D)
		dydt_TxD(2)= 0
		dydt_TxD(3)= (1+e)*tr_D
    end if
		end subroutine dydtTxD

! EULER-EXPLICIT INTEGRATION OF THE TRIAXIAL TEST
		subroutine intE_TxD
	!DESCRIPTION:
	! Perform the explicit Euler integration to simulate oedometers test using
	! dydtTxD subroutine

	HX_TxD=0
	do k=1,Number_TxDTest
		Dtime=time_TxD(k)/N_Step
		do i=1,N_POP
		! Model parameters
		phi=POP(i,1)
		hs =POP(i,2)
		n  =POP(i,3)
		ec0=POP(i,4)
		alpha=POP(i,5)
		beta=POP(i,6)
		ei0=POP(i,7)*ec0
		ed0=POP(i,8)*ec0
		! State variable
		T1=-init_TxD(k,1)/1000
		T2=-init_TxD(k,2)/1000
		e = init_TxD(k,3)

		tr_T=tr(T1,T2)
		e_d=e_dci(ed0,tr_T,hs,n)
		e_c=e_dci(ec0,tr_T,hs,n)
		e_i=e_dci(ei0,tr_T,hs,n)
		a=co_a(phi)
		f_d=fd(tr_T,e,e_d,e_c,alpha)
		f_s=fs(tr_T,e,e_d,e_c,e_i,ed0,ec0,ei0,a,hs,n,alpha,beta)
			do j=1,N_Step
			HX_TxD(j,1,i,k)=T1
			HX_TxD(j,2,i,k)=T2
			HX_TxD(j,3,i,k)=e
			! Check - tr_T < 0 e ei<e<ed
			if (tr_T > 0) exit
			if (e>e_i) exit
			if (e<e_d) exit
			call dydtTxD
			T1=dydt_TxD(1)*Dtime+T1
			T2=dydt_TxD(2)*Dtime+T2
			e =dydt_TxD(3)*Dtime+e
			tr_T=tr(T1,T2)
			e_d=e_dci(ed0,tr_T,hs,n)
			e_c=e_dci(ec0,tr_T,hs,n)
			e_i=e_dci(ei0,tr_T,hs,n)
			f_d=fd(tr_T,e,e_d,e_c,alpha)
			f_s=fs(tr_T,e,e_d,e_c,e_i,ed0,ec0,ei0,a,hs,n,alpha,beta)
			end do
		end do
	end do
		end subroutine intE_TxD

! ROSCOE VARIABLE OF A TRIAXIAL DRAINED TEST
		subroutine elab_TxD
		                   ! time_TxD,HX_TxD,RX_TxD)
	!DESCRIPTION:
	!Elaboration of the triaxial drained test respos (T1 [MPa],...

	RX_TxD=0
	do k=1,Number_TxDTest
		do i=1,N_POP
			do j=1,N_Step
			RX_TxD(j,1,i,k)=HX_TxD(j,2,i,k)-HX_TxD(j,1,i,k)
			RX_TxD(j,2,i,k)=(HX_TxD(j,3,i,k)-HX_TxD(1,3,i,k))&
			                      /(1+HX_TxD(1,3,i,k))
			end do
		end do
	end do
	do k=1,Number_TxDTest
	Dtime=time_TxD(k)/N_Step
		do i=1,N_POP
			do j=2,N_Step
			RX_TxD(j,3,i,k)=RX_TxD(j-1,3,i,k)+Dtime
			end do
		end do
	end do
		end subroutine elab_TxD

! COST FUNCTION EVALUATION IN OEDOMETRIC PLANE
		subroutine COST_EDO
	!DESCRIPTION:
	!Compute the cost function contribution in the oedometric plane

	mask_EDO_dat=.FALSE.
	do i =1,Number_EDOTest
		where (X_EDO(:,3)==real(i))
			mask_EDO_dat(:,i)=.TRUE.
		end where
		Yadm(i)=-log((1.0+MINVAL(X_EDO(:,1),mask_EDO_dat(:,i)))&
					&/(1.0+init_EDO(i,3)))
		Xadm(i)=MAXVAL(X_EDO(:,2),mask_EDO_dat(:,i))
	end do

 	do i =1,EDO_len
		k=int(X_EDO(i,3))
		Apoint_EDO(i,1)=X_EDO(i,2)/Xadm(k)
		Apoint_EDO(i,2)=-log((1+X_EDO(i,1))/(1+init_EDO(k,3)))/Yadm(k)
	end do

	do k=1,Number_EDOTest
		do i=1,N_POP
		Acurve_EDO(:,1,i,k)=-HX_EDO(:,1,i,k)/Xadm(k)
		Acurve_EDO(:,2,i,k)=-log((1+HX_EDO(:,3,i,k))/(1+init_EDO(k,3)))&
		                   &/Yadm(k)
		end do
	end do

	do k =1,Number_EDOTest
		row=MINLOC(X_EDO(:,2),mask_EDO_dat(:,k))
		mask_EDO_dat(row,k) = .FALSE.
	end do

	err_EDO=0
	do k =1,Number_EDOTest
		N_mask_Point=count(mask_EDO_dat(:,k))
		allocate(mask_Point(N_mask_Point,2))
		mask_Point(:,1)=pack(Apoint_EDO(:,1), mask_EDO_dat(:,k))
		mask_Point(:,2)=pack(Apoint_EDO(:,2), mask_EDO_dat(:,k))

		do i=1,N_POP
			err_EDO(i)=err_EDO(i)+Frechet_dist(N_Step,N_mask_Point,mask_Point,Acurve_EDO(:,:,i,k))
		end do
		deallocate(mask_Point)
	end do
		end subroutine COST_EDO

! COST FUNCTION EVALUATION IN TRIAXIAL PLANE
		subroutine COST_TxD
	!DESCRIPTION:
	!Compute the cost function contribution in the triaxials planes

	mask_TxD_dat=.FALSE.
	do i =1,Number_TxDTest
		where (X_TxD(:,4)==real(i))
			mask_TxD_dat(:,i)=.TRUE.
		end where
		Yadm_q(i)=MAXVAL(abs(X_TxD(:,3)),mask_TxD_dat(:,i))
		Xadm_q(i)=MAXVAL(abs(X_TxD(:,1)),mask_TxD_dat(:,i))
		Yadm_e(i)=MAXVAL(abs(X_TxD(:,2)),mask_TxD_dat(:,i))
		Xadm_e(i)=Xadm_q(i)
	end do
 	do i =1,TxD_len
		k=int(X_TxD(i,4))
		Apoint_TxD_q(i,1)=X_TxD(i,1)/Xadm_q(k)
		Apoint_TxD_q(i,2)=X_TxD(i,3)/Yadm_q(k)
		Apoint_TxD_e(i,1)=X_TxD(i,1)/Xadm_e(k)
		Apoint_TxD_e(i,2)=X_TxD(i,2)/Yadm_e(k)
	end do
	do k=1,Number_TxDTest
		do i=1,N_POP
		Acurve_TxD_q(:,1,i,k)=RX_TxD(:,3,i,k)*100.0/Xadm_q(k)
		Acurve_TxD_q(:,2,i,k)=RX_TxD(:,1,i,k)/Yadm_q(k)
		Acurve_TxD_e(:,1,i,k)=RX_TxD(:,3,i,k)*100.0/Xadm_e(k)
		Acurve_TxD_e(:,2,i,k)=RX_TxD(:,2,i,k)*100.0/Yadm_e(k)
		end do
	end do

	err_TxD_q=0
	do k =1,Number_TxDTest
		N_mask_Point=count(mask_TxD_dat(:,k))
		allocate(mask_Point(N_mask_Point,2))
		mask_Point(:,1)=pack(Apoint_TxD_q(:,1), mask_TxD_dat(:,k))
		mask_Point(:,2)=pack(Apoint_TxD_q(:,2), mask_TxD_dat(:,k))
		do i=1,N_POP
			err_TxD_q(i)=err_TxD_q(i)+Frechet_dist(N_Step,N_mask_Point,mask_Point,Acurve_TxD_q(:,:,i,k))
		end do
		deallocate(mask_Point)
	end do
	err_TxD_e=0

	do k =1,Number_TxDTest
		N_mask_Point=count(mask_TxD_dat(:,k))
		allocate(mask_Point(N_mask_Point,2))
		mask_Point(:,1)=pack(Apoint_TxD_e(:,1), mask_TxD_dat(:,k))
		mask_Point(:,2)=pack(Apoint_TxD_e(:,2), mask_TxD_dat(:,k))
		do i=1,N_POP
			err_TxD_e(i)=err_TxD_e(i)+Frechet_dist(N_Step,N_mask_Point,mask_Point,Acurve_TxD_e(:,:,i,k))
		end do
		deallocate(mask_Point)
	end do
		end subroutine COST_TxD

end module Evaluate_Pop
