module Testing
using DynamicalBilliards, Test
using DynamicalBilliards: isperiodic

export tag, omnibilliard, finitehorizon, testparticles, isperiodic
export acclevel

CONCAVE = Union{Semicircle}

tag(p::Particle) = "[STRAIGHT]"
tag(p::MagneticParticle) = "[MAGNETIC]"

function tag(bd::Billiard)
    s = DynamicalBilliards.isperiodic(bd) ? "[PERIODIC]" : ""
    if length(unique(map(x -> typeof(x), bd))) > 7
        s *= "[OMNI]"
    elseif any(o -> typeof(o) <: CONCAVE, bd)
        s *= "[CONCAVE]"
    end
    return s
end

tag(t::Union{<: RaySplitter, <: Tuple}) = "[RAYSPLIT]"

tag(args...) = join([tag(a) for a in args])

acclevel(::Particle{T}) where {T} = eps(T)^(3/4)
acclevel(::MagneticParticle{T}) where {T} = eps(T)^(3/4)
acclevel(bd::Billiard{T}) where {T} = isperiodic(bd) ? eps(T)^(3/5) : eps(T)^(3/4)
acclevel(args...) = maximum(acclevel(a) for a in args)

function omnibilliard()
    pref = 3 # how big is the frame with respect to the Disks
    h = pref*1.0; α = pref*0.8; r = 0.3; off = pref*0.25

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
    push!(frame, Semicircle((sp+ep)/2, α/2, normal))


    for i in 2:6
        s = iseven(i) ? β : α
        x = mod(i, 4)
        T = if x == 1
            RandomWall
        elseif x == 2
            FiniteWall
        elseif x == 3
            SplitterWall
        elseif x == 0
            InfiniteWall
        end
        sp = i != 2 ? frame[i-1].ep : startloc + α*dirs[1]
        ep = sp + s*dirs[i]
        normal = (w = ep .- sp; [-w[2], w[1]])
        push!(frame, T(sp, ep, normal))
    end


    # Radii of circles that compose the Julia logo
    offset = [0.0, off]

    R = [cos(2π/3) -sin(2π/3);
         sin(2π/3)  cos(2π/3)]


    green = Disk(center_of_mass .+ offset, r, "green")
    red = Antidot(center_of_mass .+ R*offset, r, "red")
    purple = RandomDisk(center_of_mass .+ R*R*offset, r, "purple")

    bd = Billiard(green, red, purple, frame...)
end

finitehorizon(r = 0.5) = billiard_hexagonal_sinai(0.5; setting = "periodic")

function testparticles()
    p = Particle(0.31, 0.51, 2π*rand())
    mp = MagneticParticle(0.31, 0.51, 2π*rand(), 0.5)
    return p, mp
end



function all_tests(f, args...)
    omni_tests(f, args...)
    periodic_tests(f, args...)
end
function omni_tests(f, args...)
    bd = omnibilliard()
    p, mp = testparticles()
    f(p, bd, args...)
    f(mp, bd, args...)
end
function periodic_tests(f, args...)
    bd = finitehorizon()
    p, mp = testparticles()
    f(p, bd, args...)
    f(mp, bd, args...)
end

"""
    billiards_testset(description, f, args...; caller = all_tests)
Wrap a testset around `caller(f, args...)` which times the result
and separates the test from others (`println`).
"""
function billiards_testset(d, f, args...; caller = all_tests)
    println(">> TEST: $(d)")
    t = time()
    @testset "$(d)" begin
        caller(f, args...)
    end
    println("Required time: $(round(time()-t, digits=3)) sec.")
    separator()
end
separator() = println("\n", "- "^40, "\n")

export all_tests, omni_tests, periodic_tests, billiards_testset

end
