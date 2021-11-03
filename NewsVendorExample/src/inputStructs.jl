mutable struct newsVendorModelParams
    paradigm::String
    initialInventory::Array
    orderingCost::Array
    holdingCost::Array
    backlogCost::Array
    orderingUpperBnd::Array
end

struct newsVendorUncertaintyParams
    numOfStages::Int
    dimUnc::Int
    lowerBound::Array
    upperBound::Array
    meanValue::Array
    varianceValue::Array
end

struct newsVendorLiftedParams
    brkptsVec::Array{Array,1}
    numBrkptsVec::Array
    setMatrixVec::Array{Array{Float64,2},1}  # lifted uncertainty set matrix
    setRHSVec::Array                         # lifted uncertainty set vector
    retractionMatrix::Array
    liftedMeanValue::Array
    liftedVarianceValue::Array
end

function liftedParamsNewsVendorStruct(ξParams, numBrkptsVec,brkptValuesVec)
   numStages = ξParams.numOfStages
   numUnc    = numStages-1

   W,h,R = lifted_uncertainty_set(numUnc, brkptValuesVec, numBrkptsVec, ξParams.lowerBound, ξParams.upperBound)

   ξ_liftedMean,ξ_liftedVariance = stochastic_arguments(ξParams.lowerBound,ξParams.upperBound,"uniform",brkptValuesVec, numBrkptsVec)

   liftedParams = newsVendorLiftedParams(brkptValuesVec,numBrkptsVec,W,h,R,ξ_liftedMean, ξ_liftedVariance)

   return liftedParams
end
