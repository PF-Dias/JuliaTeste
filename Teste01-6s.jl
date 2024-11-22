include("Base06.jl")        # chama o arquivo com a função que gera a base
include("Ha_Intra06.jl")    # chama o arquivo com o Hamiltoniano intra
include("Ha_Extra06.jl")    # chama o arquivo com o Hamiltoniano extra
include("Ha_Corre06.jl")    # chama o arquivo com o Hamiltoniano correção

s=zeros(Float64,6,64)
Ha_Intra=zeros(Float64,64)
Ha_Extra=zeros(Float64,64)
Ha_Corre=zeros(Float64,64)

J1=1.0
m1=1.0
T=0.00001

base(s)                    # chama a função que gera a base (6 sítios com spin de Ising)

for i in 1:64
    println(s[:,i])         # printa as configurações, linha por linha
end

print("J2=")                        # solicita o valor de J2
resposta= readline()
J2= parse(Float64, resposta)

print("J3=")                        # solicita o valor de J3
resposta= readline()
J3= parse(Float64, resposta)

#print("Passo de temperatura=")      # solicita o passo de temperatura
#resposta= readline()
#passoT= parse(Float64, resposta)

Ha_Intra6(Ha_Intra)
Ha_Extra6(Ha_Extra)
Ha_Corre6(Ha_Corre)

for i in 1:64
    println(Ha_Intra[i]+Ha_Extra[i]+Ha_Corre[i])
end

#while T ≤ 3.2
#    global T
#    Ha_Extra6(Ha_Extra)

#    Ha_Corre6(Ha_Corre)

#    println(T)
#    T=T+passoT
#end
