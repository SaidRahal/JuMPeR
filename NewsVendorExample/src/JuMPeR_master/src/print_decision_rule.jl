#-----------------------------------------------------------------------
#  print affine decision rule solution
#  Zukui Li, 2017.12.13
#-----------------------------------------------------------------------

function print_decision_rule(m::Model)

  fullDecisionRule = Dict{String, String}()
  interceptSol = Dict{String, Float64}()

  # access the robust optimization model extension
  rmext = m.ext[:JuMPeR];


  # ==============================================================
  # following code demonstrate how the adaptive variables were created in JuMPeR
  #=
  tmpM = RobustModel()

  i=1; # I arbitrarily pick one adaptive variable to demonstarte
  deps = rmext.adp_arguments[i]   # I randomly chose one for demonstration

  aux_aff = Dict{Any,Variable}()
  for j in eachindex(deps)
      # Construct a roughly sensible name for the auxiliary using
      # the adaptive variable and uncertain parameter names
      vname = string(adp_str(m,i), "{", deps[j], "}")
      aux_aff[j] = Variable(tmpM, -Inf, +Inf, :Cont, vname)
  end

  # Build the policy: intercept term
  vname = string(adp_str(m,i), "{_}")
  aux_con = Variable(tmpM, -Inf, +Inf, :Cont, vname)
  aff_policy = UncExpr(1) * aux_con

  # Build the policy: linear terms
  for j in eachindex(deps)
      push!(aff_policy, UncExpr(deps[j]), aux_aff[j])
  end
  =#
  #println("Example decision rule:")
  #println(aff_policy)

  # ======================================================================
  # here, I demonstrate how to retrive the decision rule solution and display policy
  # first, create a dictionary to store the variable name and solution values
  solutionDict = Dict(m.colNames[i] => m.colVal[i] for i in eachindex(m.colNames))
  println("Decision rule solution for adjustable variables:")
  # then, ook up the dictionary, get the value and print the affine policy
  for i in 1:length(rmext.adp_arguments)
    deps = rmext.adp_arguments[i]
    # get intercept term
    vname = string(adp_str(m,i), "{_}")
  	aux_intercept_value = solutionDict[vname]
  	# if(abs(aux_intercept_value) > 1e-6)
  		aux_intercept_value = round(aux_intercept_value;digits=3)
      	policystring = string(aux_intercept_value)
  	# else
  		# policystring = ""
  	# end

    # get linear terms
    for j in eachindex(deps)
    	vname = string(adp_str(m,i), "{", deps[j], "}")
    	aux_coeff_value = round(solutionDict[vname]; digits=3)
    	if aux_coeff_value > 0
    		policystring = policystring * " +" * string(aux_coeff_value) * "*" * string(deps[j])
    	end
    	if aux_coeff_value < 0
    		policystring = policystring * " " * string(aux_coeff_value) * "*" * string(deps[j])
    	end
    end
    fullDecisionRule[rmext.org_names[i]] = policystring
	interceptSol[rmext.org_names[i]] = aux_intercept_value
	if(policystring != "" && policystring != "0.0")
      println(rmext.org_names[i] *"= " * policystring)
    end
  end
  # return fullDecisionRule, interceptSol

end
