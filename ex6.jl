using MAT, PyPlot, LinearAlgebra
include("arnoldi.jl")
include("gramschmidt.jl")

Bwedge = matread("/Users/Jens/Documents/Documents - Jens’s MacBook Pro/PhD/Num linear algebra course/Homework 1/our code/Bwedge.mat");

B = Bwedge["B"]
B_eigvals = Bwedge["B_eigvals"]

### PART 1 (question a) and b)): plot eigenvalues and convergence approximation with disk ###

# The first eigenvalue is the one we expect fastest convergence for
ei = B_eigvals[1]
centx = -10
center = (centx,0)
centComp = centx + 0im
radius = 18

convFactor = radius/abs(centComp-ei)



scatter(real(B_eigvals),imag(B_eigvals))
plt.gcf().gca().add_artist(plt.Circle(center, radius, fill=false))
#xlim(-50,10)
#ylim(-50,10)

PyPlot.display_figs()

print("Convergence factor")
print(convFactor)


### PART 2 (question c)): error in Arnoldi method ###
m = 40
b=randn(size(B)[1]);
b = b/norm(b)
b = complex(b)

Q,H=arnoldi(B,b,m, "cgsd", false, true);

checkpoints = [2, 4, 8, 10, 20, 30, 40]

eigErrors = zeros(m)
for k=1:m

    # Arnoldi method
    Hk = H[1:k,1:k]
    F = eigen(Hk)

    if k in checkpoints
        scatter(k*ones(k), real(F.values), marker="*", color="green")
    end

    absVs = abs.(F.values .- ei)
    eigError = minimum(absVs)
    eigErrors[k] = eigError


end


PyPlot.display_figs()


semilogy(eigErrors)
PyPlot.display_figs()



### PART 3 (question d)) ###
ei = -9.8 + 2im

mu = -10
Bshift = B - mu*I

bshift = Bshift\b
bshift = bshift/norm(bshift)

Q,H=arnoldi(Bshift,bshift,m, "cgsd", true, true);

#scatter(real(eigvalues),imag(eigvalues))

#PyPlot.display_figs()


eigErrors = zeros(m)

for k=1:m

    # Arnoldi method
    Hk = H[1:k,1:k]
    F = eigen(Hk)

    eigvalues = 1 ./ F.values
    eigvalues = mu .+ eigvalues

    absVs = abs.(eigvalues .- ei)
    eigError = minimum(absVs)
    eigErrors[k] = eigError

end

plot(eigErrors)
PyPlot.display_figs()
