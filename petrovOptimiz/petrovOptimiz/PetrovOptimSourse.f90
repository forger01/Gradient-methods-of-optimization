PROGRAM test
    
    real k_4,k_5,t_3,t_4
    real t(4),tg(4),alpha
    common/min/t,tg ,alpha         
    
   omega =0.0
    !k_4=1.0
    !k_5=1.0
    !t_3=0.35
    !t_4=0.0085
    h = 2.0
    
          iter = 0
          !WRITE (6,*)a,b
          t(1) = 0.025 !k4
          t(2) = 0.0225 !k5
          t(3) = 0.35
          t(4) = 0.0085
         
     !_________________________     
         
call down()
          
          READ (5,*)QQQ

   
end

 !вычисление значения производной на каждом шаге
subroutine rp(x,y,f)
 
    
    
    real t1,t2,t3,t4,e,kds,kdm,kn,tn,k01,k02,t01,t02,y(5),f(5)
    common t1,t2,t3,t4,e,kds,kdm,kn,tn,k01,k02,t01,t02
    !tee = k4
         !write(6,*)"_____________________________", tee
    !write(6,*) f(1)
    
    f(1)=y(5)*k01/t01-y(1)/t01 !dZ
    f(2)=y(1)*(k02/t02)-2*sin(50*x)*(k02/t02)-y(2)/t02 !dY
    f(3)=10*sin(0.5*x)/t4-y(2)*(kds/t4) !dC2
    f(4)=y(3)/t2-y(1)*(kdm/t2)+10*sin(0.5*x)*t3/(t2*t4)-y(2)*(t3*kds)/(t2*t4) !dE2
    f(5)=y(4)*kn/tn-y(5)/tn+y(3)*(kn*t1)/(tn*t2)-y(1)*(kdm*kn*t1)/(tn*t2)+10*sin(0.5*x)*(t3*kn*t1)/(t4*tn*t2)-y(2)*(kds*t3*kn*t1)/(t4*tn*t2)
    
    return
end  

!модуль интегрирования методом рунге-кутта 4 порядка
    subroutine rk4(n,x,h,y,rp)
    
    real y(5),y0(5),y1(5),f(5)
    h1=0.00000
    h2=h/2
    do 11 i=1,n
    y0(i)=y(i)
11  y1(i)=y(i)
    do 12 j=1,4
    call rp(x+h1,y,f)
    h1=h2
    if(j==3)h1=h
    do 12 i=1,n
    q=h1*f(i)
    y(i)=y0(i)+q
    if(j==2) q=q+q
12  y1(i)=y1(i)+q/3.0
    do 13 i=1,n
13  y(i)=y1(i)
    !eps = 10*sin(0.5*x) - 0.1 * Y(3)
    !funk=funk+(eps*eps)*0.0001+0.1**((eps-e1)**2)/0.0001
    !e1=eps
    return
    end



!функция вычисленя функционала
real function fnkcl(k_4,k_5,t_3,t_4)
     
    external rp
    
    real t1,t2,t3,t4,e,kds,kdm,kn,tn,k01,k02,t01,t02,funk,y(5),eps,e1,k_4,k_5,t_3,t_4
    common t1,t2,t3,t4,e,kds,kdm,kn,tn,k01,k02,t01,t02
    
    !k_4=1.0
    !k_5=1.0
    !t_3=0.35
    !t_4=0.0085
    e=2.0
    kds=0.1
    kdm=1.5
    k01=1.5
    k02=0.7
    t01=0.025
    t02=0.35
    tn=0.003
    kn=25.0
    x=0.0
    x9=6.28
    h=0.00005
    kof=0.1 !коэффициент к в формуле вычисления фуекционала
    e1=0.0
    t1=k_4
    t2=k_5
    t3=t_3
    t4=t_4
    eps=0.0
!    write (6,*) k_5
    funk=0.0
    y=(/0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001/)
!1   type*,'x,x9,h,y(1),y(2),y(3),y(4)?'
    !accept*,x,x9,h,y(1),y(2),y(3),y(4)
2   call rk4(5,x,h,y,rp)
 

    x=x+h
    !open(1, file= "info.txt")
    !write(1,*)  x,y(3),funk,h1
    !close(1)
    !write(6,*),x,y(1),y(2),y(3),y(4)
    eps = 10*sin(0.5*x) - kds*y(2)
    funk=funk+(eps*eps*h)+kof*((eps-e1)**2)/h ! вычисление функционала
    e1=eps
    if((x<x9)==(h>0))goto 2
    fnkcl=funk
    !write(6,*),x,y(1),y(2),y(3),y(4)
    !write(6,*) tn
    return
    end function
   
 
   
   real function f1(alpha)
   
    
 !real alpha
          real t(4),tg(4),tee
          common/min/t,tg
         !tee = t(1)
         !write(6,*)"_____________________________", tee
          q = fnkcl(t(1)-alpha*tg(1),t(2)-alpha*tg(2),t(3)-alpha*tg(3),t(4)-alpha*tg(4))
          f1 = q
          
          return
end
   
   !минимизация функции одной переменно f1(alpha) методом половинного деления
subroutine min_1(f1, a, b)
 external f1    
    external fnkcl
   
    external grad
    external down
    external rk4
    external rp
real a1, b1, x1, x2, y1, y2,x,b2,a2
    real t(4),tg(4),alpha
    common/min/t,tg,alpha
    
    integer cnt 
    cnt =0
    a2 = a
    b2 = b
    !write (6,*)a,b
    !write (6,*)a2,b2
    
   
        
4   continue
    a1 = a2
    b1= b2
    !write(*,*)x ,a1, b1 
    eps = 0.01*(b1-a1) !задание точности вычислений
    
    do while (abs(a1-b1)>eps.and.cnt<10)
        x=(a1+b1)/2
        !write(*,*)x
            if (f1(a1)*f1(x)<0)then 
                b1=x
            else
                a1=x
            end if
cnt=cnt+1
!write(*,*)x ,a1, b1 
end do
x=(a1+b1)/2.0
  
    !x = (a1+b1)/2.0
    if((x <= a2+2.0*eps).and. (cnt < 40)) then
    a2 = a2 - 50.0*eps
    b2 = b2 - 50.0*eps
    cnt = cnt + 1
    go to 4
    end if
    if((x >= b2-2.0*eps).and. (cnt < 40)) then
    a2 = a2 + 50.0*eps
    b2 = b2 + 50.0*eps
    cnt = cnt + 1
    go to 4
    end if
    !write(6,*)a,b,x
    !min_1 = x
    !write(6,*)"x mi", x
    alpha = x
    !write(6,*)"alp mi", alpha
    return
end


!вычисление градиента   
subroutine grad(fnkcl, t_11,t_22,t_33,t_44)
REAL t_(4),dt,a,b
    real t(4),tg(4),alpha
    common/min/t,tg,alpha
    
    t_(1)=t_11
    t_(2)=t_22
    t_(3)=t_33
    t_(4)=t_44
    
    dt=0.00005
    a=fnkcl(t_(1),t_(2),t_(3),t_(4))
    do i=1,4,1
    t_(i) = t_(i)+dt
    !write (6,*) t_(i)
    b=fnkcl(t_(1),t_(2),t_(3),t_(4))
    !b=fnkcl(t_(1),t_(2),t_(3),t_(4))
    t_(i) = t_(i)-dt
    
    tg(i)= (b-a)/dt
    
    !write (6,*)"a b",a,b
     !write (6,*)"tg1 tg2 tg3 tg4",tg(1),tg(2),tg(3),tg(4)
end do
END
    
    
    
    
    
subroutine down()
     external f1    
    external fnkcl
    external min_1
    external grad
    
    external rk4
    external rp
    integer schet
    real a,b,buf,x_k(4),x_k1(4),d_k(4),d_k1(4),r_k1(4),r_k(4),r0(4),d0(4),x0(4),omega,vihod
    real t(4),tg(4),alpha
    common/min/t,tg,alpha
    schet = 0
    omega = 0.0
    vihod = 5.0
    !alpha = 0.0
     !write (6,*)"t1 t2 t3 t4",t(1),t(2),t(3),t(4)
   ! write (6,*)"tg1 tg2 tg3 tg4",tg(1),tg(2),tg(3),tg(4)
    
    !write (6,*)"fnkcl",fnkcl(t(1),t(2),t(3),t(4))
  
    
 do while((vihod>=0.01).and.(fnkcl(t(1),t(2),t(3),t(4))<10))
   call grad(fnkcl,t(1),t(2),t(3),t(4))
    WRITE (6,*)"--1--"
   WRITE (6,*)"funkcional:",fnkcl(T(1),T(2),T(3),T(4))
   WRITE (6,*)"--2--"
   WRITE (6,*)"T1:",T(1)
   WRITE (6,*)"T2:",T(2)
   WRITE (6,*)"T3:",T(3)
   WRITE (6,*)"T4:",T(4)
   
   WRITE (6,*)"--3--"
   WRITE (6,*)"grdient:",Tg(1),Tg(2),Tg(3),Tg(4)
   WRITE (6,*)"--4--"
   WRITE (6,*)"antigrdient:",-Tg(1),-Tg(2),-Tg(3),-Tg(4)
   WRITE (6,*)"--5--"
   WRITE (6,*)"A B:",a,b
   WRITE (6,*)"--6--"
   WRITE (6,*)"alpha:",alpha
   WRITE (6,*)"_______________"
   WRITE (6,*)"_______________"
   do i = 1,4
        r0(i) = tg(i)
        d0(i) = tg(i)
   end do
   
   
   do ii = 1,4
    !write (6,*) t(1)
    a=-huge(a)
    b=huge(b)
    !write (6,*) a
    
   do i = 1,4
   !write (6,*) t(i),tg(i)
        buf = abs(0.3*t(i)/tg(i))
        do while(b > buf)
            b = b-(b/2)
        end do
        do while(a < -buf)
            a = a-(a/2)
            !write (6,*) a
        end do
    end do
  end do  
   !WRITE (6,*)"funct",fnkcl(T(1),T(2),T(3),T(4)) 
  !write (6,*) a,b  
  !alpha = min_1(f1,a,b)
    call min_1(f1,a,b)
    !write (6,*) alpha
      
    do i = 1,4
        x_k(i) = t(i)-alpha*tg(i)
    end do
    
     !write (6,*) fnkcl(x_k(1),x_k(2),x_k(3),x_k(4))

   do i = 1,4
        t(i) = x_k(i)
        
        
   end do
   schet = schet+1
  
    
    vihod = tg(1)**2+tg(2)**2+tg(3)**2+tg(4)**2
    !write (6,*) vihod
   !end do
   end do
 
   
   !______________________________
   
            !WRITE (6,*)"T1 T2 T3 T4",T(1),T(2),T(3),T(4)
            !WRITE (6,*)"funct",fnkcl(T(1),T(2),T(3),T(4))
            WRITE (6,*)schet
            WRITE (6,*)"konec"

   !write (6,*)"11111111111111111111111111111111111111111111111"
    !write(6,*)t(1),t(2),fnkcl(t(1),t(2),t(3),t(4))
    !write(6,*)t(3),t(4),fnkcl(t(1),t(2),t(3),t(4))
    return
end