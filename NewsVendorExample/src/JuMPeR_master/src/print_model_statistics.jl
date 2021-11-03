#-----------------------------------------------------------------------
#  print optimization model size statistics
#  Zukui Li, 2018.1.10
#-----------------------------------------------------------------------

function print_model_statistics(mymodel::Model)

  allVarType = MathProgBase.getvartype(internalmodel(mymodel))
  numVar = MathProgBase.numvar(mymodel)
  numVarBin = count(x->x==:Bin,allVarType)
  numVarCont = count(x->x==:Cont,allVarType)

  println("Solution time: ", getsolvetime(mymodel))
  println("Objective value: ", getobjectivevalue(mymodel))
  println("Number of variables: total $numVar, binary $numVarBin, continuous $numVarCont")
  println("Number of constraints: ", MathProgBase.numconstr(mymodel))
  
end
