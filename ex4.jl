using LinearAlgebra, MatrixDepot, Random, PyPlot # for: eig, norm, etc
include("arnoldi.jl")
include("gramschmidt.jl")


print("Run! \n")


# Matrix to use for homework
nn=10;
Random.seed!(1);
A=matrixdepot("wathen",nn,nn);
n = size(A)[1]
print("Size of A\n")
print(n)
print("\n")

b=randn(size(A)[1]);
b = b/norm(b)

m=80;

Km = zeros((n,m))

# We generate Km and H for m, then we take partitions of those matrices

for k=1:m
    nb = A^(k-1)*b
    Km[:,k] = nb/norm(nb)
end

Q,H=arnoldi(A,b,m, "cgsd", false, false);



for k=1:m

    # Km method
    Kk = Km[:,1:k]
    B = Kk'*A*Kk
    C = Kk'*Kk

    F = eigen(B,C)

    # Arnoldi method
    Hk = H[1:k,1:k]
    F2 = eigen(Hk)

    print(F2.values[1])
    print("\n")

    if k%5==0
        if k==5 # just an ugly way to only print legend once
            scatter(k*ones(k), real(F.values), marker="*", color="green", label="Eigenvalue approx from (2)")
            scatter(k*ones(k), real(F2.values), marker="o", color="red", label="Eigenvalue approx from Arnoldi")
        else
            scatter(k*ones(k), real(F.values), marker="*", color="green")
            scatter(k*ones(k), real(F2.values), marker="o", color="red")
        end
    end

end


print("\n")
print("Done!")
ylim(0,500)
legend()
PyPlot.display_figs()
