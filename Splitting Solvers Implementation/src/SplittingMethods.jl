__precompile__()


module SplittingMethods


    using Reexport
    @reexport using SciMLBase
    using DiffEqBase

    using Parameters
    using OrdinaryDiffEq
    using LinearAlgebra

    include("Splitting_Coefficients.jl")
    include("Splitting_Solver.jl")
    include("Splitting_Step_Functions.jl")

    export solver_Splitting  
    export Base_3flows_mstep!
    export Base_2flows_mstep!
    export SplittingProblem
    export Splitting_alg

end # module
