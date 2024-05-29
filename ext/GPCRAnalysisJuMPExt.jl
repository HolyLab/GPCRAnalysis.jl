module GPCRAnalysisJuMPExt

using GPCRAnalysis
using JuMP
using SCS
using LinearAlgebra

function GPCRAnalysis.optimize_weights(forces)
    model = Model(SCS.Optimizer)
    set_silent(model)
    @variable(model, 0 <= weights[1:size(forces[1], 2)])
    @variable(model, t[1:length(forces)] >= 0)
    for i = eachindex(forces)
        @constraint(model, [t[i]; forces[i]*weights] in SecondOrderCone())
    end
    # @objective(model, Min, weights'*(F*weights))
    @objective(model, Min, sum(t))
    @constraint(model, sum(weights) == 1)
    optimize!(model)
    return value.(weights)
end

end
