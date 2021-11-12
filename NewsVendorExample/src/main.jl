include("startup.jl")  # if you are running main from Julia REPL, execute this line only once.

# import the parameters of the test case

solver, modelParams, uncertaintyParams = testCase()

## Affine policy
model  = newsVendorAffinePolicy(solver, modelParams, uncertaintyParams)
status = solve(model)
println("Optimal cost using affine policy is $(getobjectivevalue(model)) \n")


## Lifted Policy (default case) => same as affine policy
numOfStages  = uncertaintyParams.numOfStages
brkptSet     = [ [] for i=1:numOfStages-1    ]
numBrkptStg  = zeros(Int,numOfStages-1)
liftedParams = liftedParamsNewsVendorStruct(uncertaintyParams,numBrkptStg,brkptSet)
model  = newsVendorLiftedPolicy(solver, modelParams, uncertaintyParams, liftedParams)
status = solve(model)
println("Lifted policy with $(brkptSet[1]) breakpoint in each stage")
println("Optimal cost is $(getobjectivevalue(model)) \n")


## Lifted Policy with breakpoint at the demand mean value
numOfStages  = uncertaintyParams.numOfStages
brkptSet     = [ [ uncertaintyParams.meanValue[i] ] for i=1:numOfStages-1    ]
numBrkptStg  = ones(Int,numOfStages-1)
liftedParams = liftedParamsNewsVendorStruct(uncertaintyParams,numBrkptStg,brkptSet)
model  = newsVendorLiftedPolicy(solver, modelParams, uncertaintyParams, liftedParams)
status = solve(model)
println("Lifted policy with $(brkptSet[1]) breakpoint in each stage")
println("Optimal cost using lifted policy is $(getobjectivevalue(model)) \n")


## Lifted Policy with breakpoint at the ordering amount limit in each stage
numOfStages  = uncertaintyParams.numOfStages
brkptSet     = [ [ modelParams.orderingUpperBnd[i] ] for i=1:numOfStages-1    ]
numBrkptStg  = ones(Int,numOfStages-1)
liftedParams = liftedParamsNewsVendorStruct(uncertaintyParams,numBrkptStg,brkptSet)
model  = newsVendorLiftedPolicy(solver, modelParams, uncertaintyParams, liftedParams)
status = solve(model)
println("Lifted policy with $(brkptSet[1]) breakpoint in each stage")
println("Optimal cost using lifted policy is $(getobjectivevalue(model))")

#Comment for Said#TESTHHSHSHHS