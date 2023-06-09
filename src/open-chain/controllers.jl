using LinearAlgebra

include("polynomial.jl")
include("kinematics-dynamics.jl")

# a controller that eliminates accelaration
# just to illustrate how to work with dynamics with Mujoco
function find_torque_for_ddθ(model, data, ddθ)
    # M ddθ+ C = τ
    # sample controller for no accelaration

    nv = model.nv
    M = numpy.zeros((nv, nv))
    mujoco.mj_fullM(model, M, data.qM) # calculate M
    # later use mj_mulM to avoid sparse matrices
    M = pyconvert(Array, M) # convert to Matrix has the same effect
    C = pyconvert(Array, data.qfrc_bias)
    τ_mujoco = M * ddθ + C
    τ_pinocchio = inverse_dynamics(data.qpos, data.qvel, ddθ)

    # @show τ, prev_τ 
    # @assert maximum(abs.(τ - prev_τ)) < 1e-2
    return τ_pinocchio
end

# a controller that dampens the speed of movements
function find_torque_damper_controller(model, data)
    τ = -0.1 * data.qvel; # is damping cofficient
    τ = pyconvert(Array, τ);
    return τ;
end

# calculates spacial jacobian of end effector
function calculate_ee_jacp_jacr(model, data)
    site_id = data.site("end_effector").id
    nv = model.nv
    jacp = numpy.zeros((3, nv))
    jacr = numpy.zeros((3, nv))
    mujoco.mj_jacSite(model, data, jacp, jacr, site_id)
    jacp = pyconvert(Array, jacp)
    jacr = pyconvert(Array, jacr)
    return jacp, jacr
end

function find_torque_pd_controller(model, data, pd)
    kp = 200
    kd = 110
    ee_pos = pyconvert(Array, data.site("end_effector").xpos)
    dq = pyconvert(Array, data.qvel)
    J = calculate_ee_jacp(model, data)
    e = pd - ee_pos
    e_dot = -J *  dq
    F = kp * e + kd * e_dot
    τ = transpose(J) * F
    return τ;
end

let last_t = 0, sum = nothing
    function get_dt(t)
        cur_t = t
        ans = cur_t - last_t
        last_t = cur_t
        return ans # in seconds
    end
    global function integral(e, t)
        # numerically integrates e over time
        # returns the current value of integral
        if (sum === nothing)
            sum = similar(e)
            fill!(sum, 0)
            get_dt(t) # set up dt
        else
            sum += get_dt(t) * e
        end
        return sum
    end
end

function simple_pid_feed_forward_tracker(model, data, pd, vd)
    t = pyconvert(Float64, data.time)

    kp = 200
    kd = 110
    ki = 0.0
    ee_pos = pyconvert(Array, data.site("end_effector").xpos)
    joint_damp = 0.0
    dq = pyconvert(Array, data.qvel)
    Jp, Jr = calculate_ee_jacp_jacr(model, data)
    e = pd - ee_pos
    e_dot = -Jp *  dq
    F = kp * e + kd * e_dot + ki * integral(e, t)
    τ = transpose(Jp) * F - joint_damp * dq
    # @show F, kp, e, kd, e_dot, ki * integral(e, t)
    # τ_forward = find_torque_static_controller(model, data)
    return τ;
end


function simple_pd_feed_forward_tracker(model, data, pd, vd)
    t = pyconvert(Float64, data.time)

    kp = 20
    kd = 28
    # critical damping here?
    ee_pos = pyconvert(Array, data.site("end_effector").xpos)
    dq = pyconvert(Array, data.qvel)
    Jp, Jr = calculate_ee_jacp_jacr(model, data)
    e = pd - ee_pos
    e_dot = vd - Jp * dq
    ddx = kp * e + kd * e_dot
    ddθ = pinv(Jp) * ddx
    return find_torque_for_ddθ(model, data, ddθ)
end


function curve_feed_forward_tracker(model, data, pd, vd)
    t = pyconvert(Float64, data.time)

    ee_pos = pyconvert(Array, data.site("end_effector").xpos)
    q = pyconvert(Array, data.qpos)
    dq = pyconvert(Array, data.qvel)
    Jp, Jr = calculate_ee_jacp_jacr(model, data)

    poly = interpolate(0, ee_pos, Jp * dq, 1, pd, vd, 5)
    
    # @assert maximum(abs.(eval(poly, 0)- ee_pos)) <= 1e-6
    # @assert maximum(abs.(eval(diff(poly), 0) - (Jp * dq))) <= 1e-6
    # @assert maximum(abs.(eval(poly, 1)- pd)) <= 1e-6
    # @assert maximum(abs.(eval(diff(poly), 1) - vd)) <= 1e-6

    T = 5 # s(t) = (1-exp(-t/T))
    poly0 = eval(poly, 0)
    dpoly0 = eval(diff(poly), 0)
    ddpoly0 = eval(diff(diff(poly)), 0)
    # x(t) = poly(s(t))
    # dx(t) = dpoly(s(t)) * ds(t)
    # ddx(t) = ddpoly(s(t)) * ds(t) * ds(t) + dpoly(s(t)) * dds(t)
    ddx = (ddpoly0 - dpoly0) / (T * T)
    ddθ = pinv(Jp) * ddx

    τ = find_torque_for_ddθ(model, data, ddθ)
    push!(θ_hist, q)
    push!(dθ_hist, dq)
    push!(ddθ_hist, ddθ)
    push!(ts_hist, t)

    return τ
end

θ_hist = []
dθ_hist = []
ddθ_hist = []
ts_hist = []

function plot_controller_logs() 
    get_index(arr, index) = [x[index] for x in arr]
    # plot(ts_hist, [
    #     numeric_integrate(ts_hist, get_index(dθ_hist, 2)),
    #     get_index(θ_hist, 2)
    # ])
    # display(plot(ts_hist, [
    #     numeric_integrate(ts_hist, get_index(ddθ_hist, 2)),
    #     get_index(dθ_hist, 2)],
    #     label=["integral ddθ ordered" "dθ observed"],
    #     title="error in ordered ddθ"))
end


# set controller here
tracker_controller = (model, data, pd, vd) -> simple_pd_feed_forward_tracker(model, data, pd, vd)
# tracker_controller = (model, data, pd, vd) -> test(model, data)
