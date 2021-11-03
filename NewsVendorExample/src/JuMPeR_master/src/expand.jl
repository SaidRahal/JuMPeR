#-----------------------------------------------------------------------
# JuMPeR  --  JuMP Extension for Robust Optimization
# http://github.com/IainNZ/JuMPeR.jl
#-----------------------------------------------------------------------
# Copyright (c) 2016: Iain Dunning
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#-----------------------------------------------------------------------
# src/adaptive/expand.jl
# Adaptive robust optimization support - pre-solve expansion
#-----------------------------------------------------------------------

any_adaptive(u::UncVarExpr) = any(isa.(u.vars,Adaptive))


function expand_adaptive(rm::Model)
    # Extract robust-specifc part of model information
    rmext = get_robust(rm)::RobustModelExt
    # If no adaptive variables, bail right away
    isempty(rmext.adp_policy) && return
    # Collect the new variables or expressions that will replace
    # the adaptive variables in the current constraints
    new_vars = Any[]
    # Collect the new constraints we'll be adding to make the
    # policies make sense
    new_cons = UncConstraint[]
    # For every adaptive variable...
    for i in 1:rmext.num_adps
        pol = rmext.adp_policy[i]  # Extract the policy
        #---------------------------------------------------------------
        if pol == :Static
            # Static policy - no dependence on the uncertain parameters
            # Create a normal variable to replace this adaptive variable
            push!(new_vars, Variable(rm,
                rmext.adp_lower[i], rmext.adp_upper[i],
                rmext.adp_cat[i], rmext.adp_names[i]) )
        #---------------------------------------------------------------
        elseif pol == :Affine
            # Extract the container of uncertain parameters that this
            # adaptive variable depends on
            deps = rmext.adp_arguments[i]
            # Ensure stage-wise dependence not being used
            #if rmext.adpStage[i] != 0
            #    error("Not implemented!")
            #end
            # Create auxiliary variables - one for each uncertain parameter,
            # and one for the independent term
            aux_aff = Dict{Any,Variable}()
            for j in eachindex(deps)
                # Construct a roughly sensible name for the auxiliary using
                # the adaptive variable and uncertain parameter names
                vname = string(adp_str(rm,i), "{", deps[j], "}")
                aux_aff[j] = Variable(rm, -Inf, +Inf, :Cont, vname)
            end
            vname = string(adp_str(rm,i), "{_}")
            aux_con = Variable(rm, -Inf, +Inf, :Cont, vname)
            # Build the policy
            aff_policy = UncExpr(1) * aux_con
            for j in eachindex(deps)
                push!(aff_policy, UncExpr(deps[j]), aux_aff[j])
            end
            #println(aff_policy)
            push!(new_vars, aff_policy)
            # Add bound constraints on the policy
            if rmext.adp_lower[i] != -Inf  # Lower bound?
                push!(new_cons, UncConstraint(aff_policy, rmext.adp_lower[i], +Inf))
            end
            if rmext.adp_upper[i] != +Inf  # Lower bound?
                push!(new_cons, UncConstraint(aff_policy, -Inf, rmext.adp_upper[i]))
            end
        #---------------------------------------------------------------
        else
            error("Unknown policy type")
        end
    end

    # Create new constraints from uncertain-and-variable constraints
    # that contained an adaptive variable
    sp_objective_expanded = false
    sp_objective_flag = false
    for uncaffcon in rmext.unc_constraints
        #-----------------------------
        # if it is stochasitc programming objecitve function constraint
        # then we will handle a(ξ)x(ξ) specially using variance data
        if sp_objective_expanded == false
            sp_objective_flag = is_SP_objective(uncaffcon)
        end
        if sp_objective_flag == true
            lhs = uncaffcon.terms       # (c0+cᵀξ) + ∑ᵢ(a0ᵢ+aᵢᵀξ)xᵢ + ∑ᵢ(b0ᵢ+bᵢᵀξ)yᵢ(ξ)
            new_lhs = UncVarExpr(lhs.constant)   #  first term: (c0+cᵀξ)
            for (coeff, var) in linearterms(lhs)
                if !isa(var, Adaptive)           # second term (static xi):  (a0ᵢ+aᵢᵀξ)xᵢ
                    new_lhs += coeff * var
                else                             # Adaptive variable term (b0ᵢ+bᵢᵀξ)(y0ᵢ + yᵢᵀξ)
                    expand_adap_var = new_vars[var.id]   # yᵢ(ξ) = y0 + yᵢᵀξ
                    b0 = coeff.constant    # b0ᵢ
                    new_lhs += b0 * expand_adap_var   # part a of third term: b0ᵢ * (y0ᵢ + yᵢᵀξ)
                    for (bij, ξj) in linearterms(coeff)   # coeff=bᵢᵀξ = ∑ⱼbᵢⱼξⱼ
                        mean_ξj = rmext.unc_expectation[ξj.id]
                        var_ξj = rmext.unc_variance[ξj.id]
                        for(cj, yij) in linearterms(expand_adap_var)   # expand_adap_var = 1*y0 + ξᵀyᵢ =  1*y0ᵢ +∑ⱼξⱼyᵢⱼ
                            if cj.constant == 1.00  # cj is 1, yij is y0ᵢ
                                new_lhs += bij * ξj * yij  #   part b of third term:  bᵢᵀξ y0ᵢ      expectation is not processed here, but in apply_reform() function
                            else   # cj is ξⱼ
                                cj_id = 0
                                for (tirivalone, uncparam) in linearterms(cj)
                                    if isa(uncparam, Uncertain)
                                        cj_id = uncparam.id
                                        break
                                    end
                                end
                                mean_cj = rmext.unc_expectation[cj_id]
                                if ξj.id == cj_id
                                    E_ξjcj = var_ξj + mean_ξj*mean_ξj
                                else
                                    cov_ξjcj = 0 # In current version, assume no correlation between uncertain parameters !!!
                                    E_ξjcj = cov_ξjcj + mean_ξj*mean_cj
                                end
                                new_lhs += bij * E_ξjcj * yij  # part c of third term: bᵢᵀξ yᵢᵀξ
                            end
                        end
                    end
                end
            end
            sp_objective_expanded = true
            sp_objective_flag = false
        else
        #-----------------------------
            # regular constraint (not objective function definition)
            # zukui comment: lb <= c(ξ) + a(ξ)x + ax(ξ) <= ub
            lhs = uncaffcon.terms
            !any_adaptive(lhs) && continue
            new_lhs = UncVarExpr(lhs.constant)    # copy c(ξ) term
            for (coeff, var) in linearterms(lhs)
            new_lhs += coeff * (isa(var, Adaptive) ? new_vars[var.id] : var)
	          end		   
        end
        push!(new_cons, UncConstraint(new_lhs, uncaffcon.lb, uncaffcon.ub))
        # Remove old constraint by emptying all fields
        lhs.vars = JuMPeRVar[]
        lhs.coeffs = UncExpr[]
        lhs.constant = UncExpr()
        if uncaffcon.lb != -Inf
            uncaffcon.lb = 0
        end
        if uncaffcon.ub != +Inf
            uncaffcon.ub = 0
        end
    end

    # Create new constraints from number-and-adaptive constraints
	# zukui comment: lb <= c + a*x(ξ) <= ub										
    for varaffcon in rmext.adapt_constraints
        lhs = varaffcon.terms
        new_lhs = UncVarExpr(lhs.constant)
        for (coeff, var) in linearterms(lhs)
            new_lhs += coeff * (isa(var, Adaptive) ? new_vars[var.id] : var)
        end
        push!(new_cons, UncConstraint(new_lhs, varaffcon.lb, varaffcon.ub))
        # We don't need to do anything with the old constraints        # zukui comment: why?
    end

    # Add these constraints to the RobustModel
    map(c->JuMP.addconstraint(rm, c), new_cons)
end

#-----------------------------------------------
# zukui add 2017.Dec.15
function is_SP_objective(con)
    lhs = con.terms
    objflag = false
    for (uaff,var) in linearterms(lhs)   # uaff  --> aᵢᵀu + bᵢ,    var --> xᵢ,    uaff.constant -->  bᵢ
        if string(var) == "ϕ_stochastic_program_obj"
            objflag = true
            break
        end
    end
    return objflag
end
											