################################################################################
# Julia implementation of entropy regularised optimal transport
# explicit multiplication
#
# Annina Lieberherr, Spring 2021
################################################################################

using Pkg
using FLoops
using LinearAlgebra

debug = true

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
    # if debug
    #     # testing
    #     plan = zeros(n, m)
    #     K = zeros(n, m)
    #     left = zeros(n,n)
    #     right= zeros(n,n)
    #     for i in 1:n
    #         for j in 1:m
    #             K[i,j] = exp(-(x[i] - y[j])^2/eps)
    #             right[j,j] = vcurrent[j]
    #         end
    #         left[i,i] = ucurrent[i]
    #     end
    #     plan = left*K*right
    # end
    return cost
end

function sinkhorn_explmul(x,y,a,b,eps,numItermax=1e8,threshold=1e-08,evalStep=500)
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
            cost = eval_cost(x,y,uold,vold,eps)
        end
            
        converged = false
        if (nit % evalStep == 0)
            residual = check_convergence(ucurrent,vcurrent,x,y,eps,a)
            if (debug)
                println("residual in iteration ", nit,": ",residual)
            end
        end
        if (residual < threshold)
            cost = eval_cost(x,y,ucurrent,vcurrent,eps)
            break
        end
        nit += 1
    end
    
    return cost, residual
end
