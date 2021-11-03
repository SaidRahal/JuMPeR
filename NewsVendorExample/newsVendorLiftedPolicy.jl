function newsVendorLiftedPolicy(solver, modelParams, uncertaintyParams, liftedParams)

   ξParams   = uncertaintyParams
   numStages = ξParams.numOfStages
   numUnc    = (numStages-1)

   numUncLifted = sum(liftedParams.numBrkptsVec .+ 1 ) # single product in each stage

   numBrkptsVec = liftedParams.numBrkptsVec
   W = liftedParams.setMatrixVec
   h = liftedParams.setRHSVec
   R = liftedParams.retractionMatrix

   m = RobustModel(solver = solver)

   # --------------------------- Uncertainty information ---------------------------------
   numUncLiftedStg = liftedParams.numBrkptsVec .+ 1 # works for uni-dimensional uncertainty in stage t
   numUncLifted = sum(numUncLiftedStg)
   # Always index your uncertainty with a numeric range even if you have one element.
   @uncertain(m, ζ[n=1:numUncLifted] , expectation = liftedParams.liftedMeanValue[n], variance = liftedParams.liftedVarianceValue[n] )

   tmp_idx = cumsum(numBrkptsVec .+ 1)
   pushfirst!(tmp_idx,0)
   @constraint(m, unc_set[i=1:numUnc, j=1:numBrkptsVec[i]+2],
      sum(W[i][j,k]*ζ[tmp_idx[i]+k] for k=1:numBrkptsVec[i]+1) >= h[i][j] )

   # --------------------------- Static variables ---------------------------------
   @variable(m, x1 >= 0 )

   # ----------------------- Adaptive variables -----------------------------------
   cumUncLifted    = cumsum(numUncLiftedStg)
   if numStages > 2
      @adaptive(m, x[t=2:numStages-1] >= 0
                  , policy=Affine, depends_on=ζ[1:cumUncLifted[t-1]] )
   end
   @adaptive(m , inventory[t=2:numStages]
               , policy=Affine, depends_on=ζ[1:cumUncLifted[t-1]] )
   @adaptive(m, storage[t=2:numStages]   >= 0
               , policy=Affine, depends_on=ζ[1:cumUncLifted[t-1]] )
   @adaptive(m, backlog[t=2:numStages]   >= 0
               , policy=Affine, depends_on=ζ[1:cumUncLifted[t-1]] )
   #
   # # ------------------------------ Constraints -----------------------------------
   @constraint(m, eq1a[t=1:1], inventory[t+1] == modelParams.initialInventory[1] + x1 - (R*ζ)[t] )
   @constraint(m, eq1b[t=2:numStages-1], inventory[t+1] ==  inventory[t] + x[t] - (R*ζ)[t] )

   @constraint(m, eq2a, x1   <= modelParams.orderingUpperBnd[1] )
   @constraint(m, eq2b[t=2:numStages-1], x[t] <= modelParams.orderingUpperBnd[t] )

   @constraint(m, eq3Inventory[t=2:numStages], storage[t] >= inventory[t] )
   @constraint(m, eq3BackLog[t=2:numStages],   backlog[t] >= -inventory[t] )

   @expression(m, orderingExpr,  sum( modelParams.orderingCost[t]*x[t] for t in 2:numStages-1) + modelParams.orderingCost[1]*x1 )
   @expression(m, storageExpr, sum( modelParams.holdingCost[t]*storage[t] for t in 2:numStages) )
   @expression(m, backlogExpr, sum( modelParams.backlogCost[t]*backlog[t] for t in 2:numStages) )

   if modelParams.paradigm == "RO" # RO objective formulation
      @variable(m, obj)
      @constraint(m, obj >= orderingExpr + storageExpr + backlogExpr )
      @objective(m, Min, obj)
   else modelParams.paradigm == "SP"    # SP objective formulation
      @variable(m, ϕ_stochastic_program_obj)  # mandatory objective name for stochastic programming
      @constraint(m, ϕ_stochastic_program_obj >= orderingExpr + storageExpr + backlogExpr)
      @objective(m, Min, ϕ_stochastic_program_obj)
   end

   # ------------------------------ Solve -----------------------------------------

   return m
end
