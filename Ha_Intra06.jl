function Ha_Intra6(Ha)
    for i in 1:64
        Ha[i]=(-J1*(s[1,i]*(s[2,i]+s[6,i])+s[3,i]*(s[2,i]+s[4,i])+s[5,i]*(s[4,i]+s[6,i]))
               -J2*(s[1,i]*(s[3,i]+s[5,i])+s[2,i]*(s[4,i]+s[6,i])+s[3,i]*s[5,i]+s[4,i]*s[6,i])
               -J3*(s[1,i]*s[4,i]+s[2,i]*s[5,i]+s[3,i]*s[6,i]))

        #println(Ha[i])
    end
end