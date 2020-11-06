
function lanczos(A,b,m)
    n=length(b);
    Q=zeros(n, m+1);
    betas = zeros(m-1)
    alphas = zeros(m)

    Q[:,1]=b/norm(b);


    for k=1:m

        v=A*Q[:,k]
        alphas[k]=Q[:,k]'*v
        v=v-alphas[k]*Q[:,k]
        if k>1
            v=v-betas[k-1]*Q[:,k-1];
        end

        if k<m
            betas[k]=norm(v)
            Q[:,k+1]=v/betas[k]
        end

    end

    #H = diagm(alphas)
    H = SymTridiagonal(alphas, betas)

    return H, Q
end
