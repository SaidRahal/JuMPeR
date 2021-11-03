function newsVendorAffinePolicy(solver, modelParams, uncertaintyParams)

   ξParams   = uncertaintyParams
   numStages = ξParams.numOfStages
   numUnc    = (numStages-1) # single product (1-dimensional uncertainty)

   m = RobustModel(solver = solver)

   # Said note Aug.1.2019: Always index your uncertainty with a numeric range
   # even if you have one element. Otherwise, you will get an error
   @uncertain(m, ξParams.lowerBound[n] <= ξ[n=1:numUnc] <= ξParams.upperBound[n] , expectation = ξParams.meanValue[n], variance = ξParams.varianceValue[n])

   # --------------------------- Static variables ---------------------------------
   @variable(m, x1 >= 0 )

   # ----------------------- Adaptive variables -----------------------------------
   if numStages > 2
      @adaptive(m, x[t=2:numStages-1] >= 0 , policy=Affine, depends_on=ξ[1:(t-1)] )
   end
   @adaptive(m, inventory[t=2:numStages] , policy=Affine, depends_on=ξ[1:(t-1)] )
   @adaptive(m, storage[t=2:numStages]   >= 0, policy=Affine, depends_on=ξ[1:(t-1)] )
   @adaptive(m, backlog[t=2:numStages]   >= 0, policy=Affine, depends_on=ξ[1:(t-1)] )
   #
   # # ------------------------------ Constraints -----------------------------------
   @constraint(m, eq1a[t=1:1], inventory[t+1] == modelParams.initialInventory[1] + x1 - ξ[t] )
   @constraint(m, eq1b[t=2:numStages-1], inventory[t+1] ==  inventory[t] + x[t] - ξ[t] )

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
      @variable(m, ϕ_stochastic_program_obj)  # mandatory objective name for SP
      @constraint(m, ϕ_stochastic_program_obj >= orderingExpr + storageExpr + backlogExpr)
      @objective(m, Min, ϕ_stochastic_program_obj)
   end

   # ------------------------------ Solve -----------------------------------------

   return m
end
