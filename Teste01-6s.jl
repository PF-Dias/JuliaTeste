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

print("J2=")                        # solicita o valor de J2
J2= parse(Float64, readline())

print("J3=")                        # solicita o valor de J3
J3= parse(Float64, readline())

print("Temperatura inical=")        # solicita a temperatura inicial
T= parse(Float64, readline())

print("Tempeatura final=")          # solicita a temperatura final
Tf= parse(Float64, readline())

print("Passo de temperatura=")      # solicita o passo de temperatura
passoT= parse(Float64, readline())

arquivo="FE_J2=$(J2)_J3=$(J3)"
open(arquivo,"w") do io
println(io, "J2=$J2 \t J3=$J3")
println(io, "T \t m1 \t Z")

tempo1=time()

base(s)                 # chama a função que gera a base (6 sítios com spin de Ising)
Ha_Intra6(Ha_Intra)

while T < Tf
    global m1,ma1,T,Ha,Ha_Extra
    local Z
    Z=0.0
    ma1=m1
    erro=1.0

    while erro > 10e-9
        global m1,ma1,Ha,Ha_Extra
        local Z
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

    println("T=",T," ","m=",m1)
    println(io,"$T \t $m1 \t $Z")

    T=T+passoT
end

tempo2=time()
tempo=tempo2-tempo1
println("Tempo de execução = $(round(tempo,digits=4)) segundos")
end