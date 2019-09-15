!Module info. from home-made code writtne by Prof. Shiomi and Prof. Maruyama

module consts
  integer,parameter :: maxc=10000, numax=500
  real(8),parameter :: an=6.02217d23,                 &
       &               pl=6.626068d-34,               &
       &               kb=1.38066d-23,                &
       &               e =1.60217733d-19,             &
       &               wmc=12.d-3/an,                 &
       &               ev=1.602177d-19,               &
       &               cm=0.0299792458d0,             &
       &               mu=1.66053886d-27,             &
       &               cv=299792458d0,                &
       &               eps=1.d-16,                    &
       &               Pi=3.141592
end module consts

module root
  use consts
  integer n_root, Ntot
  integer flag(maxc), str_n(maxc)
  real mass, dt, dl0
  real root_loc(maxc,2),PE, KE 
  real Temp, Temp_r, Nl, strain(maxc),angle(maxc),strain1(maxc)
end module root

module var
  use consts
  integer Ntot,natom, time
  real*8 ru(3,maxc,numax), dr(3,maxc,numax), dv(3,maxc,numax), &
       & r(3,maxc,numax), f(3,maxc,numax), v(3,maxc,numax),    &
       & f0(3,maxc,numax), &
       & dr11(numax),dr22(numax),dr33(numax)
  real(8) k_sp, k_ang, loca,k_sp0, delt, para1, para2 
end module var


program CGMD

  use consts
  use root
  use var

  implicit integer (A-Z)

  open(111,file='../input/disp_10.avg',status='old')
  open(222,file='../input/velo_10.avg',status='old')

  open(11,file='../input/parameter.dat',status='old')
  open(12,file='../input/coord_str.dat',status='old')

  open(20,file='../result/deltT_ang_leng.dat')
  open(32,file='../result/dispM_x.dat')
  open(33,file='../result/dispM_z.dat')
  open(34,file='../result/dispM_y.dat')

  open(22,file='../result/out.dat')



!! info of the beads-spring - total number of beads
  read(12,*) natom



!! Initial location of the string, in case there are many of them
!! n_root = # of string in total 
  do i=1,n_root
	read(12,*) root_loc(i,1), root_loc(i,2)
        str_n(i)=0
  end do



!! info of a node and computation condition
  read(11,*) Nt    ! Total number of iteration 
  read(11,*) mass  ! mass of a node
  read(11,*) dt    ! time step 
  read(11,*) dl0   ! distance btw nodes
 



!! integration constant for Velocity-Verlet 
  delt=9.644*dt/mass*1000/2     



!! Spring constants and cross correlation diffusion parameter  
  k_sp0=220
  k_ang=2610/dl0

  para2=1.31                                 
  para1=2.14 


close(11)
close(12)




!! Initiation 
  do j=1,12
	call initial
	r(1,1,j)=0
	r(2,1,j)=0
	r(3,1,j)=0
	v(1,1,j)=0
	v(2,1,j)=0
	v(3,1,j)=0
  end do


!! make a bead composed of 60 carbon atom, i.e. a beads with 3 unit cell  
  do j=2,12
	do jj=1,6
		read(111,*) a1,a2,a3
		r(1,1,j)=r(1,1,j)+a1/6
		r(2,1,j)=r(2,1,j)+a2/6
		r(3,1,j)=r(3,1,j)+a3/6
		read(222,*) a1,a2,a3
		v(1,1,j)=v(1,1,j)+a1/6
		v(2,1,j)=v(2,1,j)+a2/6
		v(3,1,j)=v(3,1,j)+a3/6
	end do 

	r(2,1,j)=r(2,1,j)+36.47448+dl0   ! adjusting initial height from given in MD simulation results
	angle(j)=0
	strain(j)=0
  end do 

  ! additional node to calculate angle of first node  
  r(1,1,1)=r(1,1,2)                       
  r(2,1,1)=0
  r(3,1,1)=r(3,1,2) 


  close(111)
  close(222)


! Inital run for calculating force 
  call calcul_F


!! Nt - total number of iteration for CGMD simulation 
  do t=1,Nt

  	time=t

  	call integration


!! OUTPUT printing 
	if (mod(time,50)==0) then
		write(32,*) r(1,1,11)
		write(33,*) r(3,1,11)
		write(34,*) r(2,1,11)
	end if


	if (mod(time,10000)==0) then  
		KE=0
		do i=1,12

			KE=KE+0.5*(v(1,1,i)**2+v(2,1,i)**2+v(3,1,i)**2)*mass*mu/ev

		end do 
		write(22,*) time, KE, PE, PE+KE 
	end if 

  end do


2000 format(e15.7,'     ',e15.7)
2001 format(I5,' ',e15.7,' ',e15.7,' ',e15.7)
2002 format(a) 
2003 format(i5)
2005 format(e15.7,'         ',e15.7)


close(12)
close(20)
close(21)
close(22)
close(23)
close(24)

end program CGMD






subroutine calcul_F

  use consts
  use root
  use var

  ! implicit integer (A-Z) 
  real*8 dx, dy, dz, dl_,dr1,dr2,dr3
  real*8 dx1, dy1, dz1, dx2, dy2, dz2
  real*8 dl1, dl2, dl3,cosB,ang_B,sinB
  real*8 T, f1_ang, f2_ang, d1_1, d1_2, d1_3, d2_1, d2_2, d2_3
  real*8 td1_1, td1_2, td1_3, ABS_td
  real*8 dd1_1, dd1_2, dd1_3, dd2_1, dd2_2, dd2_3, dd3_1, dd3_2, dd3_3
  real*8 dd1_10, dd1_20, dd1_30, dd2_10, dd2_20, dd2_30, uss,vss,add_u,add_v,C
  real*8 ABS_d1, ABS_d2, L_lj, epsilon, sig, ABS_d3, test,test1
  real*8 eps_, sig_,Sum_u_s2, Sum_v_s2, Sum_uv,d_lj,add_angle2
  integer i,j

  PE=0

!
  do i=1,n_root
	do j=1,str_n(i)
 		f(1,i,j)=0
 		f(2,i,j)=0
		f(3,i,j)=0
	end do 

 	do j=2,str_n(i)-1
		dx=r(1,i,j+1)-r(1,i,j)
		dy=r(2,i,j+1)-r(2,i,j)
		dz=r(3,i,j+1)-r(3,i,j)

		dl_=sqrt(dx*dx+dy*dy+dz*dz)

		test=(dy-dl0+(dx**2+dz**2)/2/dl0)/dl0   ! Calculating Green-Lagrangian strain
		strain(j)=test
 		k_sp=k_sp0*2

 		PE=PE+test*test*k_sp*0.5*dl0	! adding PE 
 
		f(1,i,j)=f(1,i,j)+k_sp*dx/dl_*test
		f(2,i,j)=f(2,i,j)+k_sp*dy/dl_*test
		f(3,i,j)=f(3,i,j)+k_sp*dz/dl_*test

		f(1,i,j+1)=f(1,i,j+1)-k_sp*dx/dl_*test
		f(2,i,j+1)=f(2,i,j+1)-k_sp*dy/dl_*test
		f(3,i,j+1)=f(3,i,j+1)-k_sp*dz/dl_*test
	end do
  end do



  !Calculting angle                                                               
  do i=1,n_root
	do j=1,str_n(i)-2

		dx1=(r(1,i,j)-r(1,i,j+1))
		dy1=(r(2,i,j)-r(2,i,j+1))
		dz1=(r(3,i,j)-r(3,i,j+1))

		dl1=sqrt(dx1*dx1+dy1*dy1+dz1*dz1)

		dx2=(r(1,i,j+2)-r(1,i,j+1))
		dy2=(r(2,i,j+2)-r(2,i,j+1))
		dz2=(r(3,i,j+2)-r(3,i,j+1))

		dl2=sqrt(dx2*dx2+dy2*dy2+dz2*dz2)

		dx3=(r(1,i,j+2)-r(1,i,j))
		dy3=(r(2,i,j+2)-r(2,i,j))
		dz3=(r(3,i,j+2)-r(3,i,j))

		dl3=sqrt(dx3*dx3+dy3*dy3+dz3*dz3)

		cosB=(dx1*dx2+dy1*dy2+dz1*dz2)/dl1/dl2

		if (cosB >1) then 
			cosB=1
		else if (cosB<-1) then 
			cosB=-1
		end if


		! Strain for two neighbors which are composing an angle 
		test=(-dy1-dl0+(dx1**2+dz1**2)/2/dl0)/dl0	
		test1=(dy2-dl0+(dx2**2+dz2**2)/2/dl0)/dl0

		sinB=sqrt(1-cosB*cosB)
		ang_B=Pi-acos(cosB)
		T=-k_ang*ang_B


		! cross-correlation diffusion process 
		add_length=-((strain(j)-test)+(strain1(j)-test1))*para1*(-1)**(j+time)

		add_angle=(angle(j)-ang_B)*para2*(-1)**(j+time)
		add_angle2=(angle(j)-ang_B)*para2*(-1)**(j+time+1)


		angle(j)=ang_B
		strain(j)=test
		strain1(j)=test1
 
		PE=PE+ang_B*ang_B*k_ang/dl0
  
		if (sinB < 0.001) then 
			sinB=0.001
		end if


		f1_ang=-T/dl1+add_length/dl1
		f2_ang=-T/dl2+add_length/dl2

		dd1_1=dx1/dl1/sinB*cosB-dx2/dl2/sinB
		dd1_2=dy1/dl1/sinB*cosB-dy2/dl2/sinB
		dd1_3=dz1/dl1/sinB*cosB-dz2/dl2/sinB

		dd2_1=dx2/dl2/sinB*cosB-dx1/dl1/sinB
		dd2_2=dy2/dl2/sinB*cosB-dy1/dl1/sinB
		dd2_3=dz2/dl2/sinB*cosB-dz1/dl1/sinB


		f(1,i,j)=f(1,i,j)+f1_ang*dd1_1-add_angle*dx1
		f(2,i,j)=f(2,i,j)+f1_ang*dd1_2-add_angle*dy1
		f(3,i,j)=f(3,i,j)+f1_ang*dd1_3-add_angle*dz1

		f(1,i,j+1)=f(1,i,j+1)-f1_ang*dd1_1-f2_ang*dd2_1+add_angle*dx1+add_angle2*dx2
		f(2,i,j+1)=f(2,i,j+1)-f1_ang*dd1_2-f2_ang*dd2_2+add_angle*dy1+add_angle2*dy2
		f(3,i,j+1)=f(3,i,j+1)-f1_ang*dd1_3-f2_ang*dd2_3+add_angle*dz1+add_angle2*dz2

		f(1,i,j+2)=f(1,i,j+2)+f2_ang*dd2_1-add_angle2*dx2
		f(2,i,j+2)=f(2,i,j+2)+f2_ang*dd2_2-add_angle2*dy2
		f(3,i,j+2)=f(3,i,j+2)+f2_ang*dd2_3-add_angle2*dz2

	end do
  end do


end subroutine 





subroutine integration 

  use consts
  use root
  use var

  ! implicit integer (A-Z)
  real dx, dy, dz, dl_,dr1,dr2,dr3
  real dx1, dy1, dz1, dx2, dy2, dz2
  real dl1, dl2, dl3,cosB,ang_B
  real T, f1_ang, f2_ang, d1_1, d1_2, d1_3, d2_1, d2_2, d2_3
  real td1_1, td1_2, td1_3, ABS_td
  real dd1_1, dd1_2, dd1_3, dd2_1, dd2_2, dd2_3, dd3_1, dd3_2, dd3_3
  real ABS_d1, ABS_d2, L_lj, epsilon, sig, ABS_d3
  real eps_, sig_,Sum_u_s2, Sum_v_s2, Sum_uv
  integer i,j


  !Integrating velocity                                                         
  call Calcul_F

  do i=1,n_root
	do j=3,str_n(i)

		r(1,i,j)=r(1,i,j)+(dt*v(1,i,j)+f(1,i,j)*delt*dt)
		r(2,i,j)=r(2,i,j)+(dt*v(2,i,j)+f(2,i,j)*delt*dt)
		r(3,i,j)=r(3,i,j)+(dt*v(3,i,j)+f(3,i,j)*delt*dt)

		v(1,i,j)=v(1,i,j)+(f(1,i,j))*delt
		v(2,i,j)=v(2,i,j)+(f(2,i,j))*delt
		v(3,i,j)=v(3,i,j)+(f(3,i,j))*delt
	end do 
  end do 

  call Calcul_F

  do i=1,n_root
	do j=3,str_n(i)

		v(1,i,j)=v(1,i,j)+(f(1,i,j))*delt
		v(2,i,j)=v(2,i,j)+(f(2,i,j))*delt
		v(3,i,j)=v(3,i,j)+(f(3,i,j))*delt

		f0(1,i,j)=f(1,i,j)	
		f0(2,i,j)=f(2,i,j)	
		f0(3,i,j)=f(3,i,j)
	end do
  end do 

end subroutine



subroutine initial

  use consts
  use root
  use var 

  implicit integer (A-Z) 

  real :: test(2,2)
  real sum_x, sum_y, sum_z, sum_2
  real T0, c1, c2, the, phi, om 
  integer SEED,k,l

  do i=1,n_root
	Natom=Natom+1
	str_n(i)=str_n(i)+1

 	r(1,i,str_n(i))=0
 	r(2,i,str_n(i))=dl0*(str_n(i)-1)
	r(3,i,str_n(i))=0

	!dr11(str_n(i))=0
 	!dr22(str_n(i))=r(2,i,str_n(i))
 	!dr33(str_n(i))=0

 	v(1,i,str_n(i))=0
 	v(2,i,str_n(i))=0
 	v(3,i,str_n(i))=0

 	f(1,i,str_n(i))=0
 	f(2,i,str_n(i))=0
 	f(3,i,str_n(i))=0
  end do

end subroutine 




subroutine create_atom

  use consts
  use root
  use var

  implicit integer (A-Z)
  real r_eps, length, omega, c1, c2,A,B, the0,AA
  real :: test(2,2) 
  
  do i=1,n_root
	Natom=Natom+1
	str_n(i)=str_n(i)+1

	r(1,i,str_n(i))=root_loc(i,1)
	r(2,i,str_n(i))=0
	r(3,i,str_n(i))=root_loc(i,2)

	v(1,i,str_n(i))=0
	v(2,i,str_n(i))=0
	v(3,i,str_n(i))=0

	f(1,i,str_n(i))=0
	f(2,i,str_n(i))=0
	f(3,i,str_n(i))=0
  end do

  k=str_n(i-1)
  l=1
end subroutine create_atom
