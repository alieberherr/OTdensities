################################################################################
# Julia implementation of entropy regularised optimal transport
# using BLAS routines
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

function sinkhorn_mul(x,y,a,b,eps,numItermax=1e8,threshold=1e-08,tau=1e5,evalStep=500)
    n = size(x)[1]
    m = size(y)[1]


    M = zeros((n,m))
    build_cost(M,x,y,n,m)
    K = exp.(-M/eps)

    uold = reshape(ones(n)/n,(n,1))
    ucurrent = reshape(ones(n)/n,(n,1))

    # vold = reshape(ones(m)/m,(m,1))
    vcurrent = reshape(ones(m)/m,(m,1))

    nit = 0
    residual = 1.
    cost = 0.

    while (nit < numItermax)
        # copy!(uold,ucurrent)
        # copy!(vold,vcurrent)

        mul!(ucurrent,K,vcurrent)
        ucurrent .= a./ucurrent
        mul!(vcurrent,transpose(K),ucurrent)
        vcurrent .= b./vcurrent

        if (any(isnan,ucurrent) || any(isnan,vcurrent))
            println("Warning: numerical errors at iteration: ", nit,", breaking.")
            break
        end

        if (nit % evalStep == 0)
            mul!(uold,K,vcurrent) # overwrite uold
            if (debug)
                # residual = sum(abs.(a - ucurrent.*uold))
                println("residual in iteration ", nit,": ",sum(abs.(a - ucurrent.*uold)))
            end
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
