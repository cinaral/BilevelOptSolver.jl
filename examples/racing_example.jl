using BilevelOptSolver
using Random
using LinearAlgebra
include("racing.jl")
road = gen_road(seed=1234);
x0 = [1.0, 0.5, 2.0, π / 2, 0.0, -0.5, 2, π / 2];
probs = setup(;
    T=10,
    Δt=0.1,
    r=1.0,
    cost_α1=1e-3,
    cost_α2=1e-4,
    cost_β=1e-1,
    drag_coef=0.2,
    u_max_nominal=1.0,
    u_max_drafting=2.5,
    draft_box_length=5.0,
    draft_box_width=2.0,
    d=2.0,
    min_speed=-1.0,
    max_heading_offset=π / 2,
    max_heading_rate=1.0
);
solve_seq_adaptive(probs, x0, road; only_want_gnep=false, only_want_sp=false)