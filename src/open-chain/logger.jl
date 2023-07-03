using Plots
using Statistics

errors = []
times = []
xeepos = []
zeepos = []
xpdpos = []
zpdpos = []
τ_hist = []

function update_log!(model, data, t)    
    pd = pyconvert(Array, model.body("tracker-ball").pos)
    ee = pyconvert(Array, data.site("end_effector").xpos)
    τ = pyconvert(Array, data.ctrl)
    push!(times, t)
    push!(errors, pd - ee)
    push!(xeepos, ee[1])
    push!(zeepos, ee[3])
    push!(xpdpos, pd[1])
    push!(zpdpos, pd[3])
    push!(τ_hist, τ)
end

function plot_logs()
    error_norm = map(errors) do e
        sqrt(sum(e .* e))
    end
    last_sample_size = round(Int, length(error_norm) * 0.3)
    @show std(error_norm[end-last_sample_size+1:end])
    @show mean(error_norm[end-last_sample_size+1:end])

    get_index(arr, index) = [x[index] for x in arr]
    display(plot(times, 
                [get_index(τ_hist, 1), get_index(τ_hist, 2), get_index(τ_hist, 3)],
                title="control τ history"))

    px = plot(times, [xeepos, xpdpos], xlabel="time", ylabel="x pos")
    pz = plot(times, [zeepos, zpdpos], xlabel="time", ylabel="z pos")
    display(plot(px, pz, layout=(2, 1)))
    display(plot(times, error_norm, xlabel="time", ylabel="error norm"))

end