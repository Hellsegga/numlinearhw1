#Classical Gram-Schmidt
function my_hw1_cgs(Q,w,k,iter)
    # write your function
    currQ = Q[:,1:k]

    h = currQ'*w
    z = w - currQ*h

    if(iter>1)
        for i=1:(iter-1)
            g = currQ'*z
            z = z - currQ*g
        end
    end


    β = norm(z)

    return h,β,z
end

# Modified Gram-Schmidt
function my_hw1_mgs(Q,w,k)
    # write your function
    z = w
    h = zeros(k)

    for i=1:k
        h[i] = Q[:,i]'*z
        z = z - h[i]*Q[:,i]
    end

    β = norm(z)

    return h,β,z
end
