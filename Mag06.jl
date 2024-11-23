function Mag6(m1,Z,Ha)
    global m1,Z
    local v
    m1=0.0
    Z=0.0
    v=minimum(Ha)

    for i in 1:64
        m1=s[1,i]*ℯ^((v-Ha[i])/T)+m1
        Z=ℯ^((v-Ha[i])/T)+Z
    end

    m1=m1/Z
end