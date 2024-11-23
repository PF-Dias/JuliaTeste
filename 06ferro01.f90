!redde hexagonal spin 1/2
!ferromagnetico
program haxagonal6
integer:: s(6,64) 
real*8::mfe,mzz,t,m1,m2,m3,m4,m5,m6,J1,J2,J3,Ha(64),h,v,Z,F,passot,U,Uf,Ent
m1=1.d0
m2=1.d0
m3=1.d0
m4=1.d0
m5=1.d0
m6=1.d0
J1=1.d0
J2=0.d0
J3=0.2d0
T=0.00001d0
h=0.d0
passot=0.1d0

call base (s)

open(unit=26, file='TesteFE')  
write(26,*) 'Temperatura',' ','Free',' ','Internal',' ','Heat',' ','Entropia',' ','Mag'

do while (T.lt.3.5d0)
v=10.d0
    
    call auto(s,T,m1,m2,m3,m4,m5,m6,h,J1,J2,J3) 
    call ham(s,J1,J2,J3,Ha,m1,m2,m3,m4,m5,m6,h) 
    
    do i=1,64
        if (Ha(i).lt.v) then
        v=Ha(i)
        end if
    enddo

    call mag(s,T,m1,m2,m3,m4,m5,m6,Ha,v,Z,U)
    
    mfe=(m1+m2+m3+m4+m5+m6)/6.d0
    
    F=v-T*dlog(Z)
    C=(U-Uf)/passot
    if (T.eq.0.0001d0) C=0.d0
    Ent=(U-F)/T
    
        write(26,*) t,F,U,C,Ent,mfe,m1,m2,m3,m4,m5,m6
        print*, t,mfe

    T=T+passot
    Uf=U
enddo
print*, 'J3=',J3,'J2=',J2,'FE'
close(26)
end

!!!-------------------AUTO CONSISTENCIA-------------------!!!
subroutine auto(s,T,m1,m2,m3,m4,m5,m6,h,J1,J2,J3)
real*8:: T,m,erro,tol,m1,m2,m3,m4,m5,m6,Ha(64),J1,J2,J3,h,ma1,ma2,ma3,ma4,ma5,ma6,v,z,U
integer:: s(6,64),i,k
erro=1.d0
tol=10.d0**(-9.d0)
ma1=m1
ma2=m2
ma3=m3
ma4=m4
ma5=m5
ma6=m6
v=10.d0
  k=0
  do while (erro.gt.tol)
  
    call ham(s,J1,J2,J3,Ha,ma1,ma2,ma3,ma4,ma5,ma6,h)
    
    do i=1,64
        if (Ha(i).lt.v) then
        v=Ha(i)
        end if
    enddo
    
    call mag(s,T,m1,m2,m3,m4,m5,m6,Ha,v,z,U)
    
    erro1=abs(m1-ma1)
    erro2=abs(m2-ma2)
    erro3=abs(m3-ma3)
    erro4=abs(m4-ma4)
    erro5=abs(m5-ma5)
    erro6=abs(m6-ma6)
     ma1=0.1d0*m1+0.9d0*ma1
     ma2=0.1d0*m2+0.9d0*ma2
     ma3=0.1d0*m3+0.9d0*ma3
     ma4=0.1d0*m4+0.9d0*ma4
     ma5=0.1d0*m5+0.9d0*ma5
     ma6=0.1d0*m6+0.9d0*ma6
    erro=max(erro1,erro2,erro3,erro4,erro5,erro6)
    k=k+1
    
   ! if(k.gt.1000) print*, k,T,erro1,erro2,erro3,erro4,erro5,erro6
  enddo

end


!!!-------------------MAGNETIZACAO-------------------!!!
subroutine mag(s,T,m1,m2,m3,m4,m5,m6,Ha,v,z,U)
implicit none
real*8:: T,m,Z,Ha(64),m1,m2,m3,m4,m5,m6,v,U
integer:: s(6,64),j
m1=0.d0
m2=0.d0
m3=0.d0
m4=0.d0
m5=0.d0
m6=0.d0
Z=0.d0
U=0.d0

do j=1,64
    m1=m1+s(1,j)*dexp((-Ha(j)+v)/T)
    m2=m2+s(2,j)*dexp((-Ha(j)+v)/T)
    m3=m3+s(3,j)*dexp((-Ha(j)+v)/T)
    m4=m4+s(4,j)*dexp((-Ha(j)+v)/T)
    m5=m5+s(5,j)*dexp((-Ha(j)+v)/T)
    m6=m6+s(6,j)*dexp((-Ha(j)+v)/T)
    Z=Z+dexp((-Ha(j)+v)/T)
    U=U+Ha(j)*dexp((-Ha(j)+v)/T)
end do
    
    m1=m1/Z
    m2=m2/Z
    m3=m3/Z
    m4=m4/Z
    m5=m5/Z
    m6=m6/Z
    U=U/Z
end

!!!-------------------HAMILTONIANO-------------------!!!
subroutine ham(s,J1,J2,J3,Ha,m1,m2,m3,m4,m5,m6,h)
implicit none
real*8::J1,J2,J3,Ha(64),m1,m2,m3,m4,m5,m6,h,ma1,ma2,ma3,ma4,ma5,ma6
integer::i,s(6,64)
ma1=m1
ma2=m2
ma3=m3
ma4=m4
ma5=m5
ma6=m6

 do i=1,64

    Ha(i)=-j1*(s(1,i)*s(2,i)+s(2,i)*s(3,i)+s(3,i)*s(4,i)+s(4,i)*s(5,i)+s(5,i)*s(6,i)+s(6,i)*s(1,i))
    Ha(i)=Ha(i)-j2*(s(1,i)*s(3,i)+s(1,i)*s(5,i)+s(2,i)*s(4,i)+s(2,i)*s(6,i)+s(3,i)*s(5,i)+s(4,i)*s(6,i))
    Ha(i)=Ha(i)-j3*(s(1,i)*s(4,i)+s(2,i)*s(5,i)+s(3,i)*s(6,i))
    Ha(i)=Ha(i)-h*(s(1,i)+s(2,i)+s(3,i)+s(4,i)+s(5,i)+s(6,i))
    Ha(i)=Ha(i)-J1*(s(1,i)*m4+s(2,i)*ma5+s(3,i)*ma6+s(4,i)*m1+s(5,i)*ma2+s(6,i)*ma3)
    Ha(i)=Ha(i)+J1*(m1*m4+m2*ma5+m3*ma6+m4*m1+m5*ma2+m6*ma3)*0.5d0
    Ha(i)=Ha(i)-J2*(s(1,i)*(m3+ma3+m5+ma5)+s(2,i)*(m4+ma4+ma6+ma6)+s(3,i)*(m1+ma1+ma5+ma5))
    Ha(i)=Ha(i)+J2*(m1*(m3+ma3+m5+ma5)+m2*(m4+ma4+ma6+ma6)+m3*(m1+ma1+ma5+ma5))*0.5d0
    Ha(i)=Ha(i)-J2*(s(4,i)*(m2+ma2+m6+ma6)+s(5,i)*(m1+ma1+ma3+ma3)+s(6,i)*(ma2+ma2+m4+ma4))
    Ha(i)=Ha(i)+J2*(m4*(m2+ma2+m6+ma6)+m5*(m1+ma1+ma3+ma3)+m6*(ma2+ma2+m4+ma4))*0.5d0
    Ha(i)=Ha(i)-J3*(s(1,i)*(ma2+ma6)+s(2,i)*(ma1+m3)+s(3,i)*(m2+ma4)+s(4,i)*(ma3+ma5)+s(5,i)*(ma4+m6)+s(6,i)*(ma1+m5))
    Ha(i)=Ha(i)+J3*(m1*(ma2+ma6)+m2*(ma1+m3)+m3*(m2+ma4)+m4*(ma3+ma5)+m5*(ma4+m6)+m6*(ma1+m5))*0.5d0
    
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
    !print*, s(1,j),s(2,j),s(3,j),s(4,j),s(5,j),s(6,j)
end do
end do
end do
end do
end do
end do

end
