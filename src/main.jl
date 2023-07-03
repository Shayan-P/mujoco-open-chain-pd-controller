include("python-imports.jl")
include("controllers.jl")
include("logger.jl")

xml_path = joinpath(joinpath(dirname(@__DIR__), "models"), "open-chain.xml")

model = mujoco.MjModel.from_xml_path(xml_path)
data = mujoco.MjData(model)

function apply_torque!(model, data, τ)
    data.ctrl = τ;
end

function get_tracker_position(t)
    r = 1
    ω = 1.0
    pd = r .* [cos(ω * t), 0, sin(ω * t)]
    vd = r * ω .* [-sin(ω * t), 0, cos(ω * t)]
    return pd, vd
end

function apply_tracker_position!(model, data, pd)
    model.body("tracker-ball").pos = pd;
end

let visualize::Bool = true # use to turn on\off visualizer
    global function get_viewer(model, data)
        if visualize
            viewer = mujoco_viewer.MujocoViewer(model, data)
            return Dict(
                :render=> ()->viewer.render(),
                :close=> ()->viewer.close(),
                :is_alive=> ()->pyconvert(Bool, viewer.is_alive))
        else
            return Dict(
                :render=> ()->nothing,
                :close=> ()->nothing,
                :is_alive=> ()->true)
        end
    end
end
viewer = get_viewer(model, data)

mujoco.mj_forward(model, data) # necassary for initialization

for i=1:10000
    if viewer[:is_alive]()
        t = pyconvert(Float64, data.time)
        pd, vd = get_tracker_position(t)
        τ = tracker_controller(model, data, pd, vd)
        apply_torque!(model, data, τ)
        apply_tracker_position!(model, data, pd)
        update_log!(model, data, t)
        mujoco.mj_step(model, data)
        viewer[:render]()
    end
end

viewer[:close]()
plot_logs()
plot_controller_logs()
