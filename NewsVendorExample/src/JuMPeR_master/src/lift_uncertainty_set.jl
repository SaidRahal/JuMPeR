
#-----------------------------------------------------------------------
#  Lift uncertainty set for piecewise linear decision rule
#  Zukui Li, 2017.12.21
#-----------------------------------------------------------------------


# get the expecataion and variance of lifted uncertain parameters
function get_mean_var(segment_lb, segment_ub, orig_distribution, orig_lb, orig_ub)

  if orig_distribution == "uniform"
    prob_density = 1/(orig_ub-orig_lb)
    segment_mean =  prob_density * ( 0.5*(segment_ub-segment_lb)^2 +(segment_ub-segment_lb)*(orig_ub-segment_ub) )

    segment_var = 999999 # to be done, this is useful when variance is needed, e.g., a(ξ)x(ξ) in SP objective

  end

  if orig_distribution == "triangular"

    # TBD

  end

  return segment_mean, segment_var

end


#-----------------------------------------------------------------------
# uncertainty set lifting through breakpoints of each uncertain parameter
function lift_uncertainty(m::Model, ξ_lb, ξ_ub, ξ_mean, ξ_variance, breakpoints)

  numPara = length(ξ_lb)
  max_r = maximum([length(breakpoints[i])+1 for i=1:numPara])  # maximum number of segments
  ζ_mean = zeros(numPara, max_r+1)
  ζ_var = zeros(numPara, max_r+1)

  # first, define the new uncertain parameters in the lifted space
  for t=1:numPara  # for each original uncertain parameter
    r = length(breakpoints[t]) + 1  # number of segments (dimension of lifted space for this parameter)
    segment_lb = [ξ_lb[t]; breakpoints[t]]  # lower bound of segments
    segment_ub = [breakpoints[t]; ξ_ub[t]]  # upper bound of segments
    for k=1:r  # for each lifted dimension
      ζ_mean[t,k], ζ_var[t,k] =  get_mean_var(segment_lb[k], segment_ub[k], "uniform", ξ_lb[t], ξ_ub[t])
    end
    ζ_mean[t,max_r+1], ζ_var[t,max_r+1] = ξ_mean[t], ξ_variance[t]  # store original parameter mean and variance in the last element
  end


  @uncertain(m, ζ[t=1:numPara, k=1:max_r+1], expectation = ζ_mean[t, k], variance = ζ_var[t, k])  # the last element represents the original parameter


  # next, add uncertainty set constraints
  for t=1:numPara  # for each original uncertain parameter
    r = length(breakpoints[t]) + 1  # number of segments (dimension of lifted space for this parameter)
    segment_lb = [ξ_lb[t]; breakpoints[t]]  # lower bound of segments
    segment_ub = [breakpoints[t]; ξ_ub[t]]  # upper bound of segments

    # define matrix V using vertices coordinate of the convex hull
    V = [ξ_lb[t]; zeros(r-1,1)]   # first column
    tmpCol = [segment_ub[1]; zeros(r-1,1)]  # second column
    V= hcat(V, tmpCol)
    for k=2:r
      tmpCol[k] = segment_ub[k] - segment_lb[k]
      V= hcat(V, tmpCol)
    end
    V = vcat(ones(1, r+1), V)   # first row are all 1
    invV = inv(V)
    for k = 1 : r+1
      @constraint(m, invV[k,1] + sum(invV[k,j+1] * ζ[t,j] for j in 1:r) >= 0 )    # V⁻¹ * (1, ξ_1, ..., ξ_r)ᵀ >= 0
    end


  # the following part need double check, how are the equality constraint without variables handled?
# caution on the equality constraint for uncertain parameter relation !!!
# currently, this seems not working
    for j = r+1 : max_r
      @constraint(m, ζ[t,j] <= 0 )  # deactivate redundant paramters due to uneven number of breakpoints
      @constraint(m, ζ[t,j] >= 0 )
    end

    # @constraint(m, ζ[t,max_r+1] == sum(ζ[t,j] for j in 1:r) )  # relation between original parameter and the lifted parameters
    # @constraint(m, ζ[t,max_r+1] == sum(ζ[t,j] for j in 1:max_r) )  # relation between original parameter and the lifted parameters
  end

  return ζ, max_r

end
