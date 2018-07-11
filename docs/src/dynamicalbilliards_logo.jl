using DynamicalBilliards, PyPlot
using DynamicalBilliards: cross2D

# %%
h = 1.0; α = 0.8; r = 0.18; off = 0.25
# function billiards_logo(h=1.0, α=0.8, r=0.18, off=0.25)

    cos6 = cos(π/6)
    β = (h - cos6*α)/cos6
    t = α + 2β
    center_of_mass = [0.0, √3*t/6]
    startloc = [-α/2, 0.0]

    # create directions of the hexagonal 6:
    hexvert = [(cos(2π*i/6), sin(2π*i/6)) for i in 1:6]
    dirs = [SVector{2}(hexvert[i] .- hexvert[mod1(i+1, 6)]) for i in 1:6]

    frame = Obstacle{Float64}[]

    sp = startloc
    ep = startloc + α*dirs[1]
    normal = (w = ep .- sp; [-w[2], w[1]])
    push!(frame, InfiniteWall(sp, ep, normal, "frame 1"))


    for i in 2:6
        s = iseven(i) ? β : α
        T = iseven(i) ? RandomWall : InfiniteWall
        sp = frame[i-1].ep
        ep = sp + s*dirs[i]
        normal = (w = ep .- sp; [-w[2], w[1]])
        push!(frame, T(sp, ep, normal, "frame $(i)"))
    end


    # Radii of circles that compose the Julia logo
    offset = [0.0, off]

    R = [cos(2π/3) -sin(2π/3);
         sin(2π/3)  cos(2π/3)]


    green = Disk(center_of_mass .+ offset, r, "green")
    red = Antidot(center_of_mass .+ R*offset, r, "red")
    purple = RandomDisk(center_of_mass .+ R*R*offset, r, "purple")

    bd = Billiard(green, red, purple, frame...)
# end
# bd = billiards_logo(1.0, 0.8, 0.18, 0.25)


# %% Raysplitting functions for the red circle:
refraction = (φ, pflag, ω) -> pflag ? 0.5φ : 2.0φ
transmission_p = (p) -> (φ, pflag, ω) -> begin
    if pflag
        p*exp(-(φ)^2/2(π/8)^2)
    else
        abs(φ) < π/4 ? (1-p)*exp(-(φ)^2/2(π/4)^2) : 0.0
    end
end
newoantidot = ((x, bool) -> bool ? -2.0x : -0.5x)
raya = RaySplitter([2], transmission_p(0.8), refraction, newoantidot)

# %% Create and animate particle:
p = randominside(bd, 3.5)

cd()
mkpath("dynamicalbilliards")
cd("dynamicalbilliards")

plot_billiard(bd)
axis("off")

animate_evolution(p, bd, 500, raya; newfig = false, savename="logo")
