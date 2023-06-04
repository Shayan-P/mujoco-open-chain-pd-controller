using MuJoCo 
using MuJoCo.PythonCall

xml_path = joinpath(@__DIR__, "open-chain.xml")

model = mujoco.MjModel.from_xml_path(xml_path)
data = mujoco.MjData(model)
viewer = mujoco_viewer.MujocoViewer(model, data)

# a controller that eliminates accelaration
# just to illustrate how to work with dynamics with Mujoco
function find_torque_static_controller(model, data)
    # M ddθ+ C = τ
    # sample controller for no accelaration
    nv = model.nv
    ddθ = numpy.zeros((nv,))
    M = numpy.zeros((nv, nv))
    mujoco.mj_fullM(model, M, data.qM) # calculate M
    # later use mj_mulM to avoid sparse matrices
    M = pyconvert(Array, M) # convert to Matrix has the same effect
    ddθ = pyconvert(Array, ddθ)
    C = pyconvert(Array, data.qfrc_bias)
    τ = M * ddθ + C
    return τ
end

# a controller that dampens the speed of movements
function find_torque_damper_controller(model, data)
    τ = -0.1 * data.qvel; # is damping cofficient
    τ = pyconvert(Array, τ);
    return τ;
end

# calculates spacial jacobian of end effector
function calculate_ee_jacp(model, data)
    site_id = data.site("end_effector").id
    nv = model.nv
    jacp = numpy.zeros((3, nv))
    jacr = numpy.zeros((3, nv))
    mujoco.mj_jacSite(model, data, jacp, jacr, site_id)
    jacp = pyconvert(Array, jacp)
    jacr = pyconvert(Array, jacr)
    return jacp
end

function find_torque_pd_controller(model, data, pd)
    kp = 200
    kd = 40
    ee_pos = pyconvert(Array, data.site("end_effector").xpos)
    dq = pyconvert(Array, data.qvel)
    J = calculate_ee_jacp(model, data)
    e = pd - ee_pos
    e_dot = -J *  dq
    F = kp * e + kd * e_dot
    τ = transpose(J) * F
    return τ;
end

function apply_torque!(model, data, τ)
    data.ctrl = τ;
    # or assign one by one
    # actuator_names = ["torque1", "torque2", "torque3"]
    # for (i, name) in enumerate(actuator_names)
    #     data.actuator(name).ctrl[0] = τ[i]
    # end
end

controller = find_torque_pd_controller # choose controller here

for i=1:10000
    if pyconvert(Bool, viewer.is_alive)
        τ = controller(model, data, [-1, 0, 1])
        apply_torque!(model, data, τ)
        mujoco.mj_step(model, data)
        viewer.render()
    else
        break 
    end
end

viewer.close()
