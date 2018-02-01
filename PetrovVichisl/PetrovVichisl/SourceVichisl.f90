!real function fnkcl(k_4,k_5,t_3,t_4)

    external rp
    real k4,k5,t3,t4,e,kds,k0,t1,t2,kn,k01,t01,k02,t02,kdm,tn,y(5),funk,eps,e1,k_4,k_5,t_3,t_4
    common k4,k5,t3,t4,e,kds,k0,t1,t2,kn,k01,t01,k02,t02,kdm,tn 

    !k_4=0.025
    !k_5=0.0225
    !t_3=0.35
    !t_4=0.0085
    
    
    k_4=2.4999978E-02
    k_5=2.2499554E-02
    t_3=0.3500170
    t_4=1.4006269E-04
    e=2.0
    t02=0.35
    k01=1.5
    kds=0.1
    kdm=1.5
    k0=10.0
    k02=0.7
    !t1=0.025
    !t2=0.0225
    tn=0.003
    kn=25.0
    t01=0.025
    x=0.0
    x9=6.28
    h=0.00005
    kof=0.1
    e1=0.0
    t1=k_4
    t2=k_5
    t3=t_3
    t4=t_4
    funk=0.0
    eps=0.0
    y=(/0.0000001, 0.0000001, 0.0000001, 0.0000001, 0.0000001/)
!1   type*,'x,x9,h,y(1),y(2),y(3),y(4)?'
    !accept*,x,x9,h,y(1),y(2),y(3),y(4)
2   call rk4(5,x,h,y,rp)
    x=x+h
    open(1, file= "OptZ.txt") !открытие файла
    write(1,*)  y(1) !запись в файл
    !close(1)
   
    eps = 10*sin(0.5*x) - kds*y(2)
   ! write(6,*) funk
    funk=funk+(eps*eps*h)+kof*((eps-e1)**2)/h
    e1=eps
    if((x<x9)==(h>0))goto 2
    !fnkcl=funk
    close(1)
    write(6,*),funk 
    write(6,*),x,y(1),y(2),y(3),y(4)
    read(*,*)
   
    end 




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
    
subroutine rp(x,y,f)

    real k4,k5,t3,t4,e,kds,k0,t1,t2,kn,k01,t01,k02,t02,kdm,tn,y(5),f(5)
    common k4,k5,t3,t4,e,kds,k0,t1,t2,kn,k01,t01,k02,t02,kdm,tn
    !write(6,*) f(1)
    f(1)=y(5)*k01/t01-y(1)/t01 !dZ
    f(2)=y(1)*(k02/t02)-2*sin(50*x)*(k02/t02)-y(2)/t02 !dY
    f(3)=10*sin(0.5*x)/t4-y(2)*(kds/t4) !dC2
    f(4)=y(3)/t2-y(1)*(kdm/t2)+10*sin(0.5*x)*t3/(t2*t4)-y(2)*(t3*kds)/(t2*t4) !dE2
    f(5)=y(4)*kn/tn-y(5)/tn+y(3)*(kn*t1)/(tn*t2)-y(1)*(kdm*kn*t1)/(tn*t2)+10*sin(0.5*x)*(t3*kn*t1)/(t4*tn*t2)-y(2)*(kds*t3*kn*t1)/(t4*tn*t2)
    !write(6,*),k4
    return
end