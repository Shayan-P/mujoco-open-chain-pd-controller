using Plots
using Statistics

errors = []
times = []
xeepos = []
zeepos = []
xpdpos = []
zpdpos = []

function update_log!(model, data, t)
    global errors, times
    push!(times, t)
    
    pd = pyconvert(Array, model.body("tracker-ball").pos)
    ee = pyconvert(Array, data.site("end_effector").xpos)

    push!(errors, pd - ee)
    push!(xeepos, ee[1])
    push!(zeepos, ee[3])
    push!(xpdpos, pd[1])
    push!(zpdpos, pd[3])
end

function plot_logs()
    error_norm = map(errors) do e
        sqrt(sum(e .* e))
    end
    last_sample_size = round(Int, length(error_norm) * 0.3)
    @show std(error_norm[end-last_sample_size+1:end])
    @show mean(error_norm[end-last_sample_size+1:end])

    px = plot(times, [xeepos, xpdpos], xlabel="time", ylabel="x pos")
    pz = plot(times, [zeepos, zpdpos], xlabel="time", ylabel="z pos")
    plot(px, pz, layout=(2, 1))
    plot(times, error_norm, xlabel="time", ylabel="error norm")
end