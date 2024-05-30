module GPCRAnalysisJuMPExt

using GPCRAnalysis
using JuMP
using HiGHS
using LinearAlgebra

function GPCRAnalysis.optimize_weights(F::AbstractMatrix)
    model = Model(HiGHS.Optimizer)
    set_silent(model)
    @variable(model, 0 <= weights[1:size(F, 1)])
    @objective(model, Min, weights'*(F*weights))
    @constraint(model, sum(weights) == 1)
    optimize!(model)
    return value.(weights)
end

end
