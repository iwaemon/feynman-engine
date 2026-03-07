module FeynmanEngine

# Core algebra
include("algebra/operators.jl")
include("algebra/wick.jl")

# Models
include("models/lattice.jl")
include("models/hubbard.jl")

# Diagrams
include("diagrams/graph.jl")
include("diagrams/generate.jl")
include("diagrams/classify.jl")
include("diagrams/rules.jl")

# Resummation
include("resummation/channels.jl")
include("resummation/bethe_salpeter.jl")

# Numerical evaluation
include("numerical/greens_function.jl")
include("numerical/matsubara.jl")

# Visualization
include("visualization/diagram_plot.jl")

end # module
