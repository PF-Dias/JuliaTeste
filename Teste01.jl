open("J2-0.0_J3+0.2", "w") do io
for x in 1:10
y = x^2
z = x - 3
write(io, "$x\t$y\t$z\n")
end
end