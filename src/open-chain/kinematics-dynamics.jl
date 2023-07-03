let
    urdf_path = joinpath(@__DIR__, "open-chain.urdf")
    model = pin.buildModelFromUrdf(urdf_path)
    data = model.createData()
    
    model.gravity.linear = numpy.array([0, 0, -9.8]) # manually set here. urdf does not have gravity description

    global function inverse_dynamics(θ, dθ, ddθ)
        q = numpy.array(θ)
        v = numpy.array(dθ)
        a = numpy.array(ddθ)
        τ = pyconvert(Array, pin.rnea(model, data, q, v, a))
        return τ
    end
end
