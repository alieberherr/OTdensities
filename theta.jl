using ArgParse
using CSV
using DataFrames

# include("sinkhorn_mul.jl")
# include("sinkhorn_mul_logstab.jl")
include("sinkhorn_explmul_logstab.jl")

const re = r"([\+\-]?\d{1,}\.?\d{0,}[eE]?[\+\-]?\d{0,})"

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--resultsfile"
            help="filename of the result"
        "--exfile"
            help="filename of the excitations"
        "--orbitals"
            help="name of the orbital directory"
        "--dir"
            help="name of general directory"
    end
    return parse_args(s)
end

function cube_to_array(filename,vals,dV,retrn)
    natoms = 0
    xlim = 0.0
    ylim = 0.0
    zlim = 0.0
    nx = 0
    ny = 0
    nz = 0
    open(filename, "r") do datfile
        for (i,line) in enumerate(eachline(datfile))
    		if i < 3
    			continue
    		elseif i==3
    			tmp = eachmatch(re,line)
                tmp = [ match.match for match in eachmatch(re, line)]
    			natoms = parse(Int64,tmp[1])
    			xlim, ylim, zlim = parse(Float64, tmp[2]), parse(Float64, tmp[3]), parse(Float64, tmp[4])
    		elseif 4 <= i && i <=6
                tmp = [ parse(Float64,match.match) for match in eachmatch(re, line)]
    			# n = parse(Int64,tmp[1])
                n = floor(Int64,tmp[1])
    			if tmp[2] != 0
    				nx = n
    				dx = tmp[2]
                    dV[1] = tmp[2]
                end
    			if tmp[3] != 0
    				ny = n
    				dy = tmp[3]
                    dV[2] = tmp[3]
                end
    			if tmp[4] != 0
    				nz = n
    				dz = tmp[4]
                    dV[3] = tmp[4]
                end
    		elseif 6 < i && i <= 6+abs(natoms)
    			continue
    		elseif i==6+abs(natoms)+1 && natoms < 0
    			continue
    		else
                append!(vals,[ parse(Float64,match.match) for match in eachmatch(re, line)])
            end
        end
    end

    grid = zeros((nx*ny*nz,3))
	for i in 1:nx
		for j in 1:ny
			for k in 1:nz
                n = k + (j-1)*nz + (i-1)*nz*ny
                grid[n,1] = xlim + (i-1)*dV[1]
                grid[n,2] = ylim + (j-1)*dV[2]
                grid[n,3] = zlim + (k-1)*dV[3]
                # append!(grid,Array{Float64}[[xlim + (i-1)*dx, ylim + (j-1)*dy, zlim + (k-1)*dz]])
            end
        end
    end

    if (retrn)
        return grid
    end
end

function wasserstein(filenamea,filenameb,eps=.5)
    # initialise empty arrays
    vals1 = zeros(0)
    vals2 = zeros(0)
    dV = zeros(3)

    println("stats for reading of orbitals")
    @time grid=cube_to_array(filenamea,vals1,dV,true)
    @time cube_to_array(filenameb,vals2,dV,false)
    println(size(grid))
    # now that all integrals are >= 0.99, the error introduced here is small but it's needed for the Sinkhorn
    vals1 .= vals1./sqrt(sum(vals1.^2))#*prod(dV))
    vals2 .= vals2./sqrt(sum(vals2.^2))#*prod(dV))
    # println("check valsums are equal: ",sum(vals1.^2)," and ",sum(vals2.^2))
    println("stats for Sinkhorn")
    cost,res =  sinkhorn_explmul_logstab(grid,grid,vals1.^2 .+ 1e-16,vals2.^2 .+ 1e-16,eps)
    # return 1.0
    return cost
end

function main()
    print("WARNING: THIS VERSION ONLY CALCULATES THE COST PORTION AND DOES NOT CORRECT FOR SELF-CORRELATION.")

    isturbomole = true

    args = parse_commandline()
    resultsfile = args["resultsfile"]
    exfile = args["exfile"]
    orbitals = args["orbitals"]
    dir = args["dir"]

    open(resultsfile, "w") do resfile
        write(resfile,"Molecule,Functional,Excitation,Theta1,Theta2\n")
        exdata = CSV.read(exfile, DataFrame)
        for i in 1:nrow(exdata)
            failed = false
            molecule = exdata[i,"Molecule"]
            functional = exdata[i,"Functional"]
            excitation = exdata[i,"Excitation"]
            ncntrib = exdata[i,"no contr"]
            if (ncntrib == 0)
                continue
            end

            versions = 1
            if occursin("Delta",excitation)
                versions = 2
            end
            Theta = zeros(versions)

            for j in 1:versions
                for k in 1:ncntrib
                    phia = exdata[i,"occ$k"]
                    phib = exdata[i,"virt$k"]
                    c = exdata[i,"contr$k"]
                    if (isturbomole)
                        phia = replace(phia, " " => "")
                        phib = replace(phib, " " => "")

                        if occursin("e",phia) && occursin("e",phib)
                            if j==1
                                phia = phia * "1"
                                phib = phib * "1"
                            end
                            if j==2
                                phia = phia * "1"
                                phib = phib * "2"
                            end
                        elseif occursin("e",phia)
                            if j==1
                                phia = phia * "1"
                            end
                            if j==2
                                phia = phia * "2"
                            end
                        elseif occursin("e",phib)
                            if j==1
                                phib = phib * "1"
                            end
                            if j==2
                                phib = phib * "2"
                            end
                        end

                        filenamea = dir * molecule * "/" * functional * "/" * orbitals * "/" * phia * ".cub"
                        filenameb = dir * molecule * "/" * functional * "/" * orbitals * "/" * phib * ".cub"
                    else
                        filenamea = dir * molecule * "/" * molecule * "_" * functional * ".mo" * phia * ".cube"
                        filenameb = dir * molecule * "/" * molecule * "_" * functional * ".mo" * phib * ".cube"
                    end
                    println("Molecule: ",molecule,", Excitation: ", excitation,", Functional: ", functional,", Phia: ", phia,", Phib: ", phib)
                    if (isfile(filenamea) && isfile(filenameb))
                        th = wasserstein(filenamea,filenameb)
                        println("Theta: ",th," with contribution: ",c)
                        Theta[j] = Theta[j] + c*th
                    else
                        println("failed: ",molecule,", ",excitation,", ",functional,", ",phia,", ",phib)
                        failed = true
                    end

                end # ncntrib
            end # versions
            if (!failed)
                 # write row
                 for i=1:versions
                     if i==2
                         write(resfile,molecule*functional*excitation*",,$(Theta[i])\n")
                     else
                         write(resfile,molecule*functional*excitation*",$(Theta[i])\n")
                     end
                 end
            end
        end
    end
end

main()
