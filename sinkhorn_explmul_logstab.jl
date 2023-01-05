################################################################################
# Julia implementation of entropy regularised optimal transport
# log stabilised + explicit multiplication
#
# Annina Lieberherr, Spring 2021
################################################################################

using Pkg
using LinearAlgebra

function expcost_multiply(result,l,r,v,eps,alpha,beta,tol=1e-10)
    fill!(result,0.0)
    Threads.@threads for i in 1:size(result)[1]
    	for j in 1:size(v)[1]
            tmp=0.
    		for k in 1:size(r)[2]
    			tmp += .5*(l[i, k] - r[j, k])^2
            end
    		result[i] += v[j]*exp((-tmp+alpha[i]+beta[j])/eps)
        end
    end
end

function update_potentials(u,v,alpha,beta,x,y,a,b,eps)
    n = size(u)[1]
    m = size(v)[1]

    fill!(u,0.0)
    expcost_multiply(u,x,y,v,eps,alpha,beta)
    u .= a./u
    expcost_multiply(v,y,x,u,eps,beta,alpha)
    v .= b./v
end

function check_convergence(alpha,beta,u,v,x,y,eps,a)
    n = size(alpha)[1]

    tmp = zeros(n)
    expcost_multiply(tmp,x,y,v,eps,alpha,beta)
    residual = sum(abs.(a.-u.*tmp))
    return residual
end

function eval_cost(x,y,u,v,alpha,beta,eps)
    n = size(u)[1]
    m = size(v)[1]

    blcost = zeros(n)
    Threads.@threads for i in 1:n
        for j in 1:m
            costij = 0.
            for k in 1:size(x)[2]
                costij += .5*(x[i, k] - y[j, k])^2
            end
            blcost[i] += costij*u[i]*exp((alpha[i] + beta[j] - costij)/eps)*v[j]
        end
    end

    cost = sum(blcost)
    return cost
end

function sinkhorn_explmul_logstab(x,y,a,b,eps,numItermax=100000,threshold=1e-08,tau=1e5,evalStep=1000)
    n = size(x)[1]
    m = size(y)[1]

    alpha = reshape(zeros(n),(n,1))
    beta = reshape(zeros(m),(m,1))
    ucurrent = reshape(ones(n)/n,(n,1))
    vcurrent = reshape(ones(m)/m,(m,1))
    # read data from file if there is any
    if (isfile("alpha.dat") && isfile("beta.dat") && isfile("u.dat") && isfile("v.dat"))
        println("initialising from previous iteration")
        f = open("alpha.dat")
        g = open("beta.dat")
        h = open("u.dat")
        k = open("v.dat")
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

    nit = 0
    residual = 1.
    cost = 0.

    while (nit <= numItermax)
        update_potentials(ucurrent,vcurrent,alpha,beta,x,y,a,b,eps)

        if maximum(abs.(ucurrent)) > tau && maximum(abs.(vcurrent)) > tau
            alpha .+= eps*log.(ucurrent)
            beta .+= eps*log.(vcurrent)
            fill!(ucurrent,1.0)
            fill!(vcurrent,1.0)
        end

        if (nit % evalStep == 0)
            f = open("alpha.dat","w")
            g = open("u.dat","w")
            for i in 1:n
                println(f,alpha[i])
                println(g,ucurrent[i])
            end
            close(f)
            close(g)
            f = open("beta.dat","w")
            g = open("v.dat","w")
            for i in 1:m
                println(f,beta[i])
                println(g,vcurrent[i])
            end
            close(f)
            close(g)
            residual = check_convergence(alpha,beta,ucurrent,vcurrent,x,y,eps,a)
            println("residual in iteration ", nit,": ",residual)
            if (residual < threshold)
                cost = eval_cost(x,y,ucurrent,vcurrent,alpha,beta,eps)
                break
            end
        end
        nit += 1
    end

    return cost, residual
end
