include("Base06.jl")   # chama o arquivo com a função que gera a base

s=zeros(Float64,6,64)   # a matriz s é a base de configurações de spin
Ha=zeros(Float64,64)

base(s)                # chama a função que gera a base (6 sítios com spin de Ising)

for i in 1:64
    println(s[:,i])     # printa as configurações, linha por linha
end

J1=1.0

print("J2=")            # solicita o valor de J2
resposta= readline()
J2= parse(Float64, resposta)

print("J3=")            # solicita o valor de J3
resposta= readline()
J3= parse(Float64, resposta)

println("Você escolheu J2=", J2, " e J3=", J3)

for i in 1:64
    Ha[i]=-J1*(s[1,i]*(s[2,i]+s[6,i])+s[3,i]*(s[2,i]+s[4,i])+s[5,i]*(s[4,i]+s[6,i])) \ 
          -J2*(s[1,i]*(s[3,i]+s[5,i])+s[2,i]*(s[4,i]+s[6,i])+s[3,i]*s[5,i]+s[4,i]*s[6,i]) \
          -J3*(s[1,i]*s[4,i]+s[2,i]*s[5,i]+s[3,i]*s[6,i])
    println(Ha[i])
end


