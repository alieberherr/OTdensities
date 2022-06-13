################################################################################
# Julia implementation of entropy regularised optimal transport
# log stabilised + explicit multiplication
#
# Annina Lieberherr, Spring 2021
################################################################################

using Pkg
using ThreadsX
using LinearAlgebra

debug = true

function expcost_multiply(result,l,r,v,eps,alpha,beta,tol=1e-10)
    for i in 1:size(result)[1]
    	for j in 1:size(v)[1]
            tmp=0.
    		for k in 1:size(r)[2]
    			tmp += (l[i, k] - r[j, k])^2
            end
    		result[i] += v[j]*exp((-tmp+alpha[i]+beta[j])/eps)
        end
    end
end

function update_potentials(unew,vnew,ucurrent,vcurrent,alpha,beta,x,y,a,b,eps)
    n = size(ucurrent)[1]
    m = size(vcurrent)[1]
    
    tmp = zeros(n)
    
    expcost_multiply(tmp,x,y,vcurrent,eps,alpha,beta)
    for i in 1:n
        unew[i] = a[i]/tmp[i]
    end
    tmp = zeros(m)
    expcost_multiply(tmp,y,x,unew,eps,beta,alpha)
    for i in 1:n
        vnew[i] = b[i]/tmp[i]
    end
end

function check_convergence(alpha,beta,ucurrent,vcurrent,x,y,eps,a)
    n = size(alpha)[1]
    m = size(beta)[1]
    
    tmp = zeros(n)
    
    residual = 0.
    expcost_multiply(tmp,x,y,vcurrent,eps,alpha,beta)  
    for i in 1:n
        residual += abs(a[i] - ucurrent[i]*tmp[i])
    end
    return residual
end

function eval_cost(x,y,u,v,alpha,beta,eps)
    n = size(u)[1] 
    m = size(v)[1]
    
    cost = 0.
    for i in 1:n
        for j in 1:m
            costij = 0.
            for k in 1:size(x)[2]
                costij += (x[i, k] - y[j, k])^2
            end
            cost += costij*u[i]*exp((alpha[i] + beta[j] - costij)/eps)*v[j]
        end
    end
    
    return cost
end

function sinkhorn_explmul_logstab(x,y,a,b,eps,numItermax=1e8,threshold=1e-08,tau=1e5,evalStep=500)
    n = size(x)[1]
    m = size(y)[1]
    
    alpha = zeros(n)
    beta = zeros(m)
    
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
        
        update_potentials(ucurrent,vcurrent,uold,vold,alpha,beta,x,y,a,b,eps)
            
        if maximum(abs.(ucurrent)) > tau && maximum(abs.(vcurrent)) > tau || nit % evalStep == 0
            alpha = alpha + eps*log.(ucurrent)
            beta = beta + eps*log.(vcurrent)
            ucurrent = reshape(ones(n),(n,1))
            vcurrent = reshape(ones(m),(m,1))
        end
        
        if (nit % evalStep == 0)
            residual = check_convergence(alpha,beta,ucurrent,vcurrent,x,y,eps,a)
            if (debug)
                println("residual in iteration ", nit,": ",residual)
            end
            if (residual < threshold)
                cost = eval_cost(x,y,ucurrent,vcurrent,alpha,beta,eps)
                break
            end
        end
            
        nit += 1
    end
    
    return cost, residual
end
