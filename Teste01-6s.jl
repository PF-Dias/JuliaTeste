include("Base06.jl")        # chama o arquivo com a função que gera a base
include("Ha_Intra06.jl")    # chama o arquivo com o Hamiltoniano intra    
include("Ha_Extra06.jl")    # chama o arquivo com o Hamiltoniano extra  
include("Mag06.jl")  

s=zeros(Float64,6,64)
Ha_Intra=zeros(Float64,64)
Ha_Extra=zeros(Float64,64)
Ha=zeros(Float64,64)

J1=1.0
m1=1.0
T=0.00001

#open("06FE_J2-0.0_J3+0.2","w") do io

base(s)       # chama a função que gera a base (6 sítios com spin de Ising)

print("J2=")                        # solicita o valor de J2
resposta= readline()
J2= parse(Float64, resposta)

print("J3=")                        # solicita o valor de J3
resposta= readline()
J3= parse(Float64, resposta)

print("Passo de temperatura=")      # solicita o passo de temperatura
resposta= readline()
passoT= parse(Float64, resposta)

Ha_Intra6(Ha_Intra)

while T ≤ 3.5
    global m1,ma1,T,Ha,Ha_Extra
    local Z
    Z=0.0
    ma1=m1
    erro=1.0

    while erro > 10e-9
        global m1,ma1,Ha,Ha_Extra
        local Z,v
        Z=0.0
        cont=0
        n=10

        Ha_Extra6(Ha_Extra,ma1)

        Ha=Ha_Intra+Ha_Extra

        Mag6(m1,Z,Ha)

        erro=abs(m1-ma1)
        ma1=m1
    end

    Ha_Extra=zeros(Float64,64)
    Ha_Extra6(Ha_Extra,m1)

    Ha=Ha_Intra+Ha_Extra

    Mag6(m1,Z,Ha)

    T=round(T,digits=5)
    m1=round(m1,digits=5)

    println("T=",T," ","m=",m1," ","Z=",Z)

    #write(io, "$T\t$m1\t$Z\n")

    T=T+passoT
end
#end