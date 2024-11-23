function Ha_Extra6(Ha,m1)
    for i in 1:64
        Ha[i]=(-J1*(s[1,i]+s[2,i]+s[3,i]+s[4,i]+s[5,i]+s[6,i])*m1 
              -J2*(s[1,i]+s[2,i]+s[3,i]+s[4,i]+s[5,i]+s[6,i])*4*m1
              -J3*(s[1,i]+s[2,i]+s[3,i]+s[4,i]+s[5,i]+s[6,i])*2*m1
              +(J1*3+J2*12+J3*6)*m1*m1)
    end
end