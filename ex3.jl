using LinearAlgebra, MatrixDepot, Random # for: eig, norm, etc
include("arnoldi.jl")
include("gramschmidt.jl")


print("Run! \n")

# First test matrix
#A=randn(100,100);

# Matrix to use for homework
nn=100;
Random.seed!(1);
A=matrixdepot("wathen",nn,nn);

print("Size of A\n")
print(size(A)[1])
print("\n")

b=randn(size(A)[1]);
#b = b/norm(b)

m=5;
# "cgs", "mgs", "cgsd", "cgst"
@time Q,H=arnoldi(A,b,m, "cgst", false, false);
#should_be_zero1=norm(Q*H-A*Q[:,1:m])
should_be_zero2=norm(Q'*Q-I)


print(should_be_zero2)
print("\n")
print("Done!")
