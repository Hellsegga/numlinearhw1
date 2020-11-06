"""
    Q,H=arnoldi(A,b,m)

A simple implementation of the Arnoldi method.
The algorithm will return an Arnoldi "factorization":
Q*H[1:m+1,1:m]-A*Q[:,1:m]=0
where Q is an orthogonal basis of the Krylov subspace
and H a Hessenberg matrix.

The function `my_hw1_gs(Q,w,k)` needs to be available.

Example:
```julia-repl
using Random
A=randn(100,100); b=randn(100);
m=10;
Q,H=arnoldi(A,b,m);
println("1:st should be zero = ", norm(Q*H-A*Q[:,1:m]));
println("2:nd Should be zero = ", norm(Q'*Q-I));
```

"""
function arnoldi(A,b,m,method, shift, complex)

    n=length(b);
    if complex==true
        Q=zeros(ComplexF64, n, m+1);
        H=zeros(ComplexF64, m+1, m);
    else
        Q=zeros(n,m+1);
        H=zeros(m+1,m);
    end

    Q[:,1]=b/norm(b);

    for k=1:m
        if shift
            w = A\Q[:,k]
        else
            w=A*Q[:,k]; # Matrix-vector product with last element
        end

        # Orthogonalize w against columns of Q.
        # Implement this function or replace call with code for orthogonalizatio
        if method=="cgs"
            h,β,z=my_hw1_cgs(Q,w,k,1);
        elseif method=="cgsd"
            h,β,z=my_hw1_cgs(Q,w,k,2);
        elseif method=="cgst"
            h,β,z=my_hw1_cgs(Q,w,k,3);
        elseif method=="mgs"
            h,β,z=my_hw1_mgs(Q,w,k)
        end
        #Put Gram-Schmidt coefficients into H
        H[1:(k+1),k]=[h;β];
        # normalize
        Q[:,k+1]=z/β;
    end
    return Q,H
end
