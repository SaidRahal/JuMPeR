# function lifted_uncertainty_set(brkpt_vec::Vector{Float64},n_brkpt::Int,
#             θ_lb::Number, θ_ub::Number)
function lifted_uncertainty_set(numOrigUncertainty, brkpt_vec,n_brkpt,θ_lb, θ_ub)

    dim_unc   = sum(n_brkpt .+1)
    W  = [ zeros(Float64,n_brkpt[i]+2,n_brkpt[i]+1) for i=1:numOrigUncertainty ]
    h  = [ zeros(Float64,n_brkpt[i]+2) for i=1:numOrigUncertainty ]
    R  =   zeros(Float64, numOrigUncertainty, dim_unc) # R'ξ′ = ξ
    for i=1:numOrigUncertainty
        if n_brkpt[i] > 0
            tmp_v1 = vcat(θ_lb[i], brkpt_vec[i][:])
            tmp_v2 = vcat(brkpt_vec[i][:], θ_ub[i])
            tmp    = 1 ./ (tmp_v2 .- tmp_v1)
            for ii=1:n_brkpt[i]+1
                W[i][(ii-1)+1:(ii-1)+2,ii] = [-tmp[ii],tmp[ii]]
            end
            h[i][1] = -brkpt_vec[i][1]*tmp[1]
            h[i][2] = θ_lb[i]*tmp[1]
        else
            W[i][1:2,1] = [-1,1]
            h[i][1:2,1] = [-θ_ub[i],θ_lb[i]]
        end
    end

    for i=1:numOrigUncertainty
        if i==1
            idx0 = 0
        else
            idx0 = sum(n_brkpt[1:i-1] .+ 1)
        end
        idxf = sum(n_brkpt[1:i] .+1)
        R[i,idx0+1:idxf] .= 1
    end

    return W,h,R
end

# get the expecataion and variance of lifted uncertain parameters
function stochastic_arguments(unc_lb, unc_ub, unc_distribution, brkpt_vecs, n_brkpt )

  num_unc             = length(unc_lb)
  unc_lifted_mean     = Vector{Float64}()
  unc_lifted_variance = Vector{Float64}()
  for i = 1:num_unc
    if unc_distribution == "uniform" # same distribution for all uncertainties
      if n_brkpt[i] == 0
        push!(unc_lifted_mean,0.5*(unc_ub[i] + unc_lb[i]) )
        push!(unc_lifted_variance, 1/(12*(unc_ub[i] -unc_lb[i])^2 ) )
      else
        segment_ub   = vcat(brkpt_vecs[i],unc_ub[i])  # will be updated for each ξ
        segment_lb   = vcat(unc_lb[i],brkpt_vecs[i])
        prob_density = 1/(unc_ub[i]-unc_lb[i])
        # specific formula for "first" segment
        tmp =  prob_density * ( 0.5*(segment_ub[1]^2-unc_lb[i]^2) + segment_ub[1]*(unc_ub[i]-segment_ub[1]) )
        push!(unc_lifted_mean,tmp)
        for k= 2:n_brkpt[i]+1
          tmp =  prob_density * ( 0.5*(segment_ub[k]-segment_lb[k])^2
               +(segment_ub[k]-segment_lb[k])*(unc_ub[i]-segment_ub[k]) )
          push!(unc_lifted_mean,tmp)
        end

        for k=1:n_brkpt[i]+1
          push!(unc_lifted_variance, 999)
        end
      end
    end
  end
  if unc_distribution == "triangular"

    # To Be Done

  end

  return unc_lifted_mean, unc_lifted_variance
end
