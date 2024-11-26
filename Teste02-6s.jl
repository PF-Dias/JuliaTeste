

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

open("FE_J2=$(J2)_J3=$(J3)","w") do io
println(io, "J2=$J2 \t J3=$J3")
println(io, "T \t m1 \t Z")

function Ha_Extra6(Ha,m1)
    for i in 1:64
        Ha[i]=(-J1*(s[1,i]+s[2,i]+s[3,i]+s[4,i]+s[5,i]+s[6,i])*m1 
              -J2*(s[1,i]+s[2,i]+s[3,i]+s[4,i]+s[5,i]+s[6,i])*4*m1
              -J3*(s[1,i]+s[2,i]+s[3,i]+s[4,i]+s[5,i]+s[6,i])*2*m1
              +(J1*3+J2*12+J3*6)*m1*m1)
    end
end

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

tempo1=time()

j=1
while j ≤ 64
    for a in -1:2:1
    for b in -1:2:1
    for c in -1:2:1
    for d in -1:2:1
    for e in -1:2:1
    for f in -1:2:1
        s[1,j]=a
        s[2,j]=b
        s[3,j]=c
        s[4,j]=d
        s[5,j]=e
        s[6,j]=f
        j=j+1
    end
    end
    end
    end
    end
    end
end

for i in 1:64
    Ha_Intra[i]=(-J1*(s[1,i]*(s[2,i]+s[6,i])+s[3,i]*(s[2,i]+s[4,i])+s[5,i]*(s[4,i]+s[6,i]))
           -J2*(s[1,i]*(s[3,i]+s[5,i])+s[2,i]*(s[4,i]+s[6,i])+s[3,i]*s[5,i]+s[4,i]*s[6,i])
           -J3*(s[1,i]*s[4,i]+s[2,i]*s[5,i]+s[3,i]*s[6,i]))
end

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