function  testCase()
   paradigm    = "SP"
   numOfStages = 4

   # parameters for a single product
   initialInventory = [4]
   orderingCost     = 3*ones(numOfStages-1)
   storageCost      = 1.5*ones(numOfStages-1)
   backlogCost      = 7*ones(numOfStages-1)
   pushfirst!(storageCost,0) # cost parameter not defined at t=1
   pushfirst!(backlogCost,0)
   orderingUpperBnd =  8*ones(numOfStages-1)
   modelParams = newsVendorModelParams(paradigm,initialInventory,orderingCost,storageCost, backlogCost,orderingUpperBnd)

   # uniform uncertainty distribution
   dimUnc        = 1 # per stage
   lowerBound    = zeros(numOfStages-1)
   upperBound    = 10*ones(numOfStages-1)
   meanValue     = 0.5*(lowerBound + upperBound)
   varianceValue = ((upperBound .- lowerBound).^2)./12

   uncertaintyParams = newsVendorUncertaintyParams(numOfStages,dimUnc,lowerBound,upperBound,meanValue,varianceValue)

   solver = GLPKSolverLP()

   return solver, modelParams, uncertaintyParams
end
