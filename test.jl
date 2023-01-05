include("sinkhorn_mul.jl")
include("sinkhorn_mul_logstab.jl")
# include("theta.jl")
include("sinkhorn_explmul.jl")
include("sinkhorn_explmul_logstab.jl")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function run(n)
    x = reshape([1:1.:n;],(n,1))
    a = reshape([1:1.:n;]/sum(x),(n,1))
    b = reshape([n:-1.:1;]/sum(x),(n,1))
    eps=0.5

    # final cost: 0.9902036908126535
    # filenamea = "test/polyacetylene2/pbe/orbs_adjgrid_spacing0_8/1bg.cub"
    # filenameb = "test/polyacetylene2/pbe/orbs_adjgrid_spacing0_8/2au.cub"
    # final cost: ?
    # filenamea = "../molecules_turbomole/tripeptide/pbe/orbs_adjgrid_spacing0_8/40a'.cub"
    # filenameb = "../molecules_turbomole/tripeptide/pbe/orbs_adjgrid_spacing0_8/11a\".cub"

    # vals1 = zeros(0)
    # vals2 = zeros(0)
    # dV = zeros(3)
    #
    # println("stats for reading of orbitals")
    # @time grid=cube_to_array(filenamea,vals1,dV,true)
    # @time cube_to_array(filenameb,vals2,dV,false)
    #
    # vals1 .= vals1./sqrt(sum(vals1.^2))#*prod(dV))
    # vals2 .= vals2./sqrt(sum(vals2.^2))#*prod(dV))
    # println(size(vals1))


    # println("explicit for-loop multiplication")
    # @time cost, residual = sinkhorn_explmul(x,x,a,b,eps)
    # println("Final cost: ", cost)
    #
    println("explicit for-loop multiplication, log-stabilised")
    @time cost, residual = sinkhorn_explmul_logstab(x,x,a,b,eps,"test")
    println("Final cost: ", cost)
    #
    # println("multiplication from Julia:")
    # @time cost, residual = sinkhorn_mul(x,x,a,b,eps)
    # @time cost, residual = sinkhorn_mul(grid,grid,vals1.^2 .+ 1e-16,vals2.^2 .+ 1e-16,eps)
    # println("Final cost: ", cost)

    # println("multiplication from Julia, log-stabilised:")
    # @time cost, residual = sinkhorn_mul_logstab(x,x,a,b,eps)
    # @time cost, residual = sinkhorn_mul_logstab(grid,grid,vals1.^2 .+ 1e-16,vals2.^2 .+ 1e-16,eps)
    # println("Final cost: ", cost)
end

run(100)
