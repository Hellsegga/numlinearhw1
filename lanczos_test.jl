using LinearAlgebra, MatrixDepot, Random, Arpack # for: eig, norm, etc
include("lanczos.jl")

print("Run! \n")

# First test matrix
#A=randn(100,100);

# Matrix to use for homework
nn=100;
Random.seed!(1);
A=matrixdepot("wathen",nn,nn);


m = 20
b=randn(size(A)[1]);
b = b/norm(b)


H = lanczos(A, b, m)


print("Done!")

#FF = eigs(A, nev=20, v0=b)
