################################################################################
# Julia implementation of entropy regularised optimal transport
# explicit multiplication
#
# Annina Lieberherr, Spring 2021
################################################################################

using Pkg
using LinearAlgebra

function expcost_multiply(result,l,r,u,eps,tol=1e-10)
    for i in 1:size(result)[1]
    	for j in 1:size(u)[1]
            tmp=0.
    		for k in 1:size(r)[2]
    			tmp += (l[i, k] - r[j, k])^2
            end
    		result[i] += u[j]*exp(-tmp/eps)
        end
    end
end

function costexpcost_multiply(result,l,r,u,eps,tol=1e-10)
    for i in 1:size(result)[1]
    	for j in 1:size(u)[1]
            tmp=0.
    		for k in 1:size(r)[2]
    			tmp += (l[i, k] - r[j, k])^2
            end
    		result[i] += u[j]*tmp*exp(-tmp/eps)
        end
    end
end

function update_potentials(unew,vnew,ucurrent,vcurrent,x,y,a,b,eps)
    n = size(ucurrent)[1]
    m = size(vcurrent)[1]

    tmp = zeros(n)

    expcost_multiply(tmp,x,y,vcurrent,eps)
    for i in 1:n
        unew[i] = a[i]/tmp[i]
    end
    tmp = zeros(m)
    expcost_multiply(tmp,y,x,unew,eps)
    for i in 1:n
        vnew[i] = b[i]/tmp[i]
    end
end

function check_convergence(ucurrent,vcurrent,x,y,eps,a)
    n = size(ucurrent)[1]
    m = size(vcurrent)[1]

    tmp = zeros(n)

    residual = 0.
    expcost_multiply(tmp,x,y,vcurrent,eps)
    for i in 1:n
        residual += abs(a[i] - ucurrent[i]*tmp[i])
    end
    return residual
end

function eval_cost(x,y,ucurrent,vcurrent,eps)
    n = size(ucurrent)[1]
    m = size(vcurrent)[1]
    MKdotv = zeros(n)
    cost = 0.
    costexpcost_multiply(MKdotv,x,y,vcurrent,eps)
    for i in 1:size(MKdotv)[1]
        cost += ucurrent[i]*MKdotv[i]
    end
    return cost
end

function sinkhorn_explmul(x,y,a,b,eps,numItermax=1e8,threshold=1e-08,evalStep=50)
    n = size(x)[1]
    m = size(y)[1]

    uold = reshape(ones(n)/n,(n,1))
    ucurrent = reshape(ones(n)/n,(n,1))

    vold = reshape(ones(m)/m,(m,1))
    vcurrent = reshape(ones(m)/m,(m,1))

    nit = 0
    residual = 1.
    cost = 0.

    while (nit <= numItermax)
        uold = copy(ucurrent)
        vold = copy(vcurrent)

        update_potentials(ucurrent,vcurrent,uold,vold,x,y,a,b,eps)

        if (any(isnan,ucurrent) || any(isnan,vcurrent))
            println("Warning: numerical errors at iteration: ", nit,", breaking", nit)
            break
        end

        converged = false
        if (nit % evalStep == 0)
            residual = check_convergence(ucurrent,vcurrent,x,y,eps,a)
            println("residual in iteration ", nit,": ",residual)
        end
        if (residual < threshold)
            cost = eval_cost(x,y,ucurrent,vcurrent,eps)
            break
        end
        nit += 1
    end

    return cost, residual
end

function sinkhorn_explmul_logstab(x,y,a,b,eps,name,numItermax=1000,threshold=1e-08,tau=1e5,evalStep=10)
    n = size(x)[1]
    m = size(y)[1]

    alpha = reshape(zeros(n),(n,1))
    beta = reshape(zeros(m),(m,1))
    # read data from file if there is any
    println("name:",name)
    ucurrent = reshape(ones(n)/n,(n,1))
    vcurrent = reshape(ones(m)/m,(m,1))
    if (isfile("alpha"*name*".dat") && isfile("beta"*name*".dat") && isfile("u"*name*".dat") && isfile("v"*name*".dat"))
        println("initialising from previous iteration")
        f = open("alpha"*name*".dat")
        g = open("beta"*name*".dat")
        h = open("u"*name*".dat")
        k = open("v"*name*".dat")
        for i in 1:n
            alpha[i] = parse(Float64,readline(f))
            ucurrent[i] = parse(Float64,readline(h))
        end
        for i in 1:m
            beta[i] = parse(Float64,readline(g))
            vcurrent[i] = parse(Float64,readline(k))
        end
        close(f)
        close(g)
        close(h)
        close(k)
    end
    ucurrent = reshape(ones(n)/n,(n,1))
    vcurrent = reshape(ones(m)/m,(m,1))

    nit = 0
    residual = 1.
    cost = 0.

    while (nit <= numItermax)
        update_potentials(ucurrent,vcurrent,alpha,beta,x,y,a,b,eps)

        if (maximum(abs.(ucurrent)) > tau || maximum(abs.(vcurrent)) > tau) || nit==numItermax
            alpha .+= eps*log.(ucurrent)
            beta .+= eps*log.(vcurrent)
            fill!(ucurrent,1.0)
            fill!(vcurrent,1.0)
        end

        if (nit % evalStep == 0)
            f = open("alpha"*name*".dat","w")
            g = open("u"*name*".dat","w")
            for i in 1:n
                println(f,alpha[i])
                println(g,ucurrent[i])
            end
            close(f)
            close(g)
            f = open("beta"*name*".dat","w")
            g = open("v"*name*".dat","w")
            for i in 1:m
                println(f,beta[i])
                println(g,vcurrent[i])
            end
            close(f)
            close(g)
            residual = check_convergence(alpha,beta,ucurrent,vcurrent,x,y,eps,a)
            println("residual in iteration ", nit,": ",residual)
            if (residual < threshold)
                alpha .+= eps*log.(ucurrent)
                beta .+= eps*log.(vcurrent)
                f = open("alpha"*name*".dat","w")
                for i in 1:n
                    println(f,alpha[i])
                end
                close(f)
                f = open("beta"*name*".dat","w")
                for i in 1:m
                    println(f,beta[i])
                end
                cost = eval_cost(x,y,ucurrent,vcurrent,alpha,beta,eps)
                break
            end
        end
        nit += 1
    end

    return cost, residual
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
