!redde hexagonal spin 1/2
!ferromagnetico
program haxagonal6
integer:: s(6,64),taxa,t1,t2
real*8::mfe,mzz,t,m1,J1,J2,J3,Ha0(64),Ha1(64),Ha(64),h,v,Z,F,passot,tempo
m1=1.d0
J1=1.d0
J2=0.d0
J3=0.2d0
T=0.0001d0
passot=0.0001d0

    ! Obtém a taxa do relógio (ticks por segundo)
    call system_clock(count_rate=taxa)

    ! Tempo inicial
    call system_clock(count=t1)

call base (s)

open(unit=26, file='TesteFE')  
write(26,*) 'Temperatura',' ','Free',' ','Internal',' ','Heat',' ','Entropia',' ','Mag'

    call ham0(s,J1,J2,J3,Ha0)

do while (T.lt.3.5d0)
v=10.d0
    
    call auto(s,T,m1,J1,J2,J3,Ha0) 
    call ham(s,J1,J2,J3,Ha1,m1) 
    
    do i=1,64
        Ha(i)=Ha0(i)+Ha1(i)
    enddo

    call mag(s,T,m1,Ha)
    
        write(26,*) t,m1
        print*, T,m1
        
    T=T+passot
    Uf=U
enddo
print*, 'J3=',J3,'J2=',J2,'FE'
close(26)

    ! Tempo final
    call system_clock(count=t2)

    ! Calcula o tempo decorrido em segundos
    tempo=abs(real(t2 - t1)/real(taxa))
    print*, "tempo total (segundos)=", tempo
end

!!!-------------------AUTO CONSISTENCIA-------------------!!!
subroutine auto(s,T,m1,J1,J2,J3,Ha0)
real*8:: T,erro,tol,m1,Ha0(64),Ha1(64),Ha(64),J1,J2,J3,ma1
integer:: s(6,64),i,k
erro=1.d0
tol=10.d0**(-9.d0)
ma1=m1

  do while (erro.gt.tol)
  
    call ham(s,J1,J2,J3,Ha1,ma1)
    
    do i=1,64
        Ha(i)=Ha0(i)+Ha1(i)
    enddo
    
    call mag(s,T,m1,Ha)
    
    erro=abs(m1-ma1)
    ma1=m1

  enddo

end


!!!-------------------MAGNETIZACAO-------------------!!!
subroutine mag(s,T,m1,Ha)
implicit none
real*8:: T,Z,Ha(64),m1,v
integer:: s(6,64),j
m1=0.d0
Z=0.d0
v=minval(Ha)

do j=1,64
    m1=m1+s(1,j)*dexp((-Ha(j)+v)/T)
    Z=Z+dexp((-Ha(j)+v)/T)
end do
    
    m1=m1/Z
end

!!!-------------------HAMILTONIANO-------------------!!!
subroutine ham0(s,J1,J2,J3,Ha)
implicit none
real*8::J1,J2,J3,Ha(64)
integer::i,s(6,64)

 do i=1,64

    Ha(i)=-J1*(s(1,i)*(s(2,i)+s(6,i))+s(3,i)*(s(2,i)+s(4,i))+s(5,i)*(s(4,i)+s(6,i)))
    Ha(i)=Ha(i)-J2*(s(1,i)*(s(3,i)+s(5,i))+s(2,i)*(s(4,i)+s(6,i))+s(3,i)*s(5,i)+s(4,i)*s(6,i))
    Ha(i)=Ha(i)-J3*(s(1,i)*s(4,i)+s(2,i)*s(5,i)+s(3,i)*s(6,i))
    
end do
 
end 

!!!-------------------HAMILTONIANO-------------------!!!
subroutine ham(s,J1,J2,J3,Ha,m1)
implicit none
real*8::J1,J2,J3,Ha(64),m1
integer::i,s(6,64)

 do i=1,64

    Ha(i)=-J1*(s(1,i)+s(2,i)+s(3,i)+s(4,i)+s(5,i)+s(6,i))*m1
    Ha(i)=Ha(i)-J2*(s(1,i)+s(2,i)+s(3,i)+s(4,i)+s(5,i)+s(6,i))*4*m1
    Ha(i)=Ha(i)-J3*(s(1,i)+s(2,i)+s(3,i)+s(4,i)+s(5,i)+s(6,i))*2*m1
    Ha(i)=Ha(i)+(J1*3+J2*12+J3*6)*m1*m1
    
end do
 
end 

!!!-------------------BASE-------------------!!!
subroutine base(s)
implicit none
integer:: s(6,64),a,b,c,d,e,f,j
j=0

do a=-1,1,2
do b=-1,1,2
do c=-1,1,2
do d=-1,1,2
do e=-1,1,2
do f=-1,1,2
    j=j+1
    s(1,j)=a
    s(2,j)=b
    s(3,j)=c
    s(4,j)=d
    s(5,j)=e
    s(6,j)=f
end do
end do
end do
end do
end do
end do

end
