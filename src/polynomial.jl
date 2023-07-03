struct Polynomial 
    cof::AbstractArray
end

function diff(p::Polynomial)
    new_cof = [p.cof[i+1] * i for i in 1:(length(p.cof)-1)]
    return Polynomial(new_cof)
end

function eval(p::Polynomial, x)
    return sum(p.cof[i+1] * (x^i) for i in 0:(length(p.cof)-1))
end

Base.string(p::Polynomial) = begin
    if length(p.cof) == 0
        return "0"
    else 
        return join(
            ["($(p.cof[i+1])*x$i)" for i in 0:(length(p.cof)-1)]
            ,"+")
    end
end

Base.print(io::IO, p::Polynomial) = print(io, Base.string(p))
Base.show(io::IO, p::Polynomial) = println(io, Base.string(p))
Base.repr(p::Polynomial) = Base.string(p)

function interpolate(x0, f0, df0, x1, f1, df1, n::Int)::Polynomial
    """an n degree polynomial that satisfies below equations:
        # P(x0) = f0
        # dP(x0) = df0
        # P(x1) = f1
        # dP(x1) = df1
    """
    A = [
        reshape([x0^i for i in 0:n], (1, n+1));
        reshape([(i == 0 ? 0 : i * (x0^(i-1))) for i in 0:n], (1, n+1));
        reshape([x1^i for i in 0:n], (1, n+1));
        reshape([(i == 0 ? 0 : i * (x1^(i-1))) for i in 0:n], (1, n+1))
    ]
    b = [f0, df0, f1, df1]
    cof = pinv(A) * b
    return Polynomial(cof)
end


function numeric_integrate(t, y)
    sm = 0
    tbef = 0
    res = []
    for (t0, y0) in zip(t, y)
        sm += (t0-tbef) * y0
        tbef = t0
        push!(res, sm)
    end
    return res
end
