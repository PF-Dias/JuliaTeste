function Auto6(m1)
    include("Ha_Extra06.jl")
    global m1
    tol=10e-9
    erro=1.0
    ma1=m1

    while erro > tol
        global m1,ma1,erro
        local Ha
        
        Ha_Extra6(Ha_Extra)

        Ha=Ha_Intra+Ha_Extra
        x=minimum(Ha)

        # call Mag

        erro=abs(m1-ma1)
        ma1=m1

    end

end