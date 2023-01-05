################################################################################
# Julia implementation of entropy regularised optimal transport
# log-stabilised + using BLAS routines
#
# Annina Lieberherr, Spring 2021
################################################################################

using Pkg
using LinearAlgebra

function build_cost(M,x,y,n,m)
    for i in 1:n
        for j in 1:m
            for k in 1:size(x)[2]
                M[i,j] += (x[i,k] - y[j,k])*(x[i,k] - y[j,k])
            end
        end
    end
end

function check_convergence(ucurrent,vcurrent,x,y,eps,a)
    n = size(ucurrent)[1]
    m = size(vcurrent)[1]
    tmp = zeros(n)

    residual = 0.
    tmp = u.*(K*vcurrent)
    for i in 1:n
        residual += abs(a[i] - tmp[i])
    end
    residual
end

function sinkhorn_mul_logstab(x,y,a,b,eps,numItermax=5,threshold=1e-08,tau=1e5,evalStep=500)
    n = size(x)[1]
    m = size(y)[1]

    M = zeros((n,m))
    build_cost(M,x,y,n,m)
    K = exp.(-M/eps)

    uold = reshape(ones(n)/n,(n,1))
    ucurrent = reshape(ones(n)/n,(n,1))

    vold = reshape(ones(m)/m,(m,1))
    vcurrent = reshape(ones(m)/m,(m,1))

    alpha = reshape(zeros(n)/n,(n,1))
    beta = reshape(zeros(n)/n,(n,1))

    nit = 0
    residual = 1.
    cost = 0.

    while (nit <= numItermax)
        copy!(uold,ucurrent)
        copy!(vold,vcurrent)

        # println("ucurrent tmp:",ucurrent)
        mul!(ucurrent,K,vcurrent)
        # println("ucurrent tmp:",ucurrent)
        ucurrent .= a./ucurrent
        mul!(vcurrent,K',ucurrent)
        vcurrent .= b./vcurrent

        # println("ucurrent:",ucurrent)

        if maximum(abs.(ucurrent)) > tau && maximum(abs.(vcurrent)) > tau || nit % evalStep == 0
            alpha .= alpha .+ eps*log.(ucurrent)
            beta .= beta .+ eps*log.(vcurrent)
            fill!(K,0.0)
            @inbounds @fastmath for i in eachindex(alpha), j in eachindex(beta) # update K
                K[i,j] = exp((-M[i,j]+alpha[i] + beta[j])/eps)
            end
            ucurrent .= reshape(ones(n),(n,1))
            vcurrent .= reshape(ones(m),(m,1))
            # println("exp(alpha/eps):",exp.(alpha/eps))
        end

        if (nit % evalStep == 0)
            mul!(uold,K,vcurrent) # overwrite uold
            println("residual in iteration ", nit,": ",sum(abs.(a - ucurrent.*uold)))
            if (sum(abs.(a - ucurrent.*uold)) < threshold)
                mul!(uold,(M.*K),vcurrent) # overwrite uold
                cost = dot(ucurrent,uold)
                break
            end
        end
        nit += 1
    end

    cost, residual
end
