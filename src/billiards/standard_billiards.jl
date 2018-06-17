using StaticArrays

export billiard_rectangle, billiard_sinai, billiard_polygon, billiard_lorentz,
billiard_raysplitting_showcase, billiard_hexagonal_sinai, billiard_bunimovich,
billiard_stadium, billiard_mushroom

####################################################
## Famous/Standard Billiards
####################################################
"""
```julia
billiard_rectangle(x=1.0, y=1.0; setting = "standard")
```
Return a vector of obstacles that defines a rectangle billiard of size (`x`, `y`).

### Settings
* "standard" : Specular reflection occurs during collision.
* "periodic" : The walls are `PeriodicWall` type,
  enforcing periodicity at the boundaries
* "random" : The velocity is randomized upon collision.
* "ray-splitting" : All obstacles in the billiard allow for ray-splitting.
"""
function billiard_rectangle(x=1.0, y=1.0; setting::String = "standard")

    x = convert(AbstractFloat, x)
    x, y = promote(x,y)
    o = typeof(x)(0.0)
    if setting == "standard"
        sp = [o,y]; ep = [o, o]; n = [x,o]
        leftw = InfiniteWall(sp, ep, n, "Left wall")
        sp = [x,o]; ep = [x, y]; n = [-x,o]
        rightw = InfiniteWall(sp, ep, n, "Right wall")
        sp = [x,y]; ep = [o, y]; n = [o,-y]
        topw = InfiniteWall(sp, ep, n, "Top wall")
        sp = [o,o]; ep = [x, o]; n = [o,y]
        botw = InfiniteWall(sp, ep, n, "Bottom wall")
    elseif setting == "periodic"
        sp = [o,y]; ep = [o, o]; n = [x,o]
        leftw = PeriodicWall(sp, ep, n, "Left periodic boundary")
        sp = [x,o]; ep = [x, y]; n = [-x,o]
        rightw = PeriodicWall(sp, ep, n, "Right periodic boundary")
        sp = [x,y]; ep = [o, y]; n = [o,-y]
        topw = PeriodicWall(sp, ep, n, "Top periodic boundary")
        sp = [o,o]; ep = [x, o]; n = [o,y]
        botw = PeriodicWall(sp, ep, n, "Bottom periodic boundary")
    elseif setting == "random"
        sp = [o,y]; ep = [o, o]; n = [x,o]
        leftw = RandomWall(sp, ep, n, "Left random wall")
        sp = [x,o]; ep = [x, y]; n = [-x,o]
        rightw = RandomWall(sp, ep, n, "Right random wall")
        sp = [x,y]; ep = [o, y]; n = [o,-y]
        topw = RandomWall(sp, ep, n, "Top random wall")
        sp = [o,o]; ep = [x, o]; n = [o,y]
        botw = RandomWall(sp, ep, n, "Bottom random wall")
    elseif setting == "ray-splitting"
        sp = [o,y]; ep = [o, o]; n = [x,o]
        leftw = SplitterWall(sp, ep, n, "Left ray-splitting wall")
        sp = [x,o]; ep = [x, y]; n = [-x,o]
        rightw = SplitterWall(sp, ep, n, "Right ray-splitting wall")
        sp = [x,y]; ep = [o, y]; n = [o,-y]
        topw = SplitterWall(sp, ep, n, "Top ray-splitting wall")
        sp = [o,o]; ep = [x, o]; n = [o,y]
        botw = SplitterWall(sp, ep, n, "Bottom ray-splitting wall")
    else
        throw(ArgumentError("The given setting=$setting is unknown."))
    end
    return Billiard(botw, rightw, topw, leftw)
end

"""
```julia
billiard_sinai(r=0.25, x=1.0, y=1.0; setting = "standard")
```
Return a vector of obstacles that defines a Sinai billiard of size (`x`, `y`) with
a disk in its center, of radius `r`.

In the periodic case, the system is also known as "Lorentz Gas".

### Settings
* "standard" : Specular reflection occurs during collision.
* "periodic" : The walls are `PeriodicWall` type,
  enforcing periodicity at the boundaries
* "random" : The velocity is randomized upon collision.
* "ray-splitting" : All obstacles in the billiard allow for ray-splitting.
"""
function billiard_sinai(r=0.25, x=1.0, y=1.0; setting = "standard")
    if (setting == "periodic") && (r>=x/2 || r>=y/2)
        es = "Disk radius too big for a periodic Sinai billiard.\n"
        es*= "Obstacles must not overlap with `PeriodicWall`s."
        error(es)
    end
    r, x, y = promote(r,x,y)
    bdr = billiard_rectangle(x, y; setting = setting)

    c = [x/2, y/2]
    if setting == "random"
        centerdisk = RandomDisk(c, r, "Random disk")
    elseif setting == "ray-splitting"
        centerdisk = Antidot(c, r, "Antidot")
    else
        centerdisk = Disk(c, r, "Disk")
    end

    return Billiard(centerdisk, bdr...)
end

"""
    billiard_lorentz(r=0.25, x=1.0, y=1.0)
Alias for `billiard_sinai(r,x,y; setting = "periodic")`.
"""
billiard_lorentz(r=0.25, x=1.0, y=1.0) = billiard_sinai(r,x,y; setting = "periodic")

"""
    billiard_polygon(n::Int, R, center = [0,0]; setting = "standard")
Return a vector of obstacles that defines a regular-polygonal billiard
with `n` sides, radius `r` and given `center`.

Note: `R` denotes the so-called outer radius, not the inner one.

### Settings
* "standard" : Specular reflection occurs during collision.
* "periodic" : The walls are `PeriodicWall` type, enforcing periodicity
  at the boundaries. Only available for `n=4` or `n=6`.
* "random" : The velocity is randomized upon collision.
"""
function billiard_polygon(sides::Int, r::Real, center = [0,0]; setting = "standard")
    S = typeof(convert(AbstractFloat, r))
    bd = Obstacle{S}[]
    verteces = [S[r*cos(2π*i/sides), r*sin(2π*i/sides)] .+ center for i in 1:sides]

    if setting == "standard"
        T = InfiniteWall
        wallname = "wall"
    elseif setting == "periodic"
        if sides != 4 && sides != 6
            error("Polygonal and periodic billiard can exist only for `n=4` or `n=6`")
        end
        T = PeriodicWall
        wallname = "periodic wall"
        inr = S(r*cos(π/sides))
    elseif setting == "random"
        T = RandomWall
        wallname = "random wall"
    end

    for i in eachindex(verteces)
        starting = verteces[i]
        ending = verteces[mod1(i+1, sides)]
        # Normal vector must look at where the particle is coming from
        w = ending - starting
        if setting == "periodic"
            normal = 2inr*normalize([-w[2], w[1]])
            wall = T(starting, ending, normal, wallname*" $i")
        else
            normal = [-w[2], w[1]]
            wall = T(starting, ending, normal, wallname*" $i")
        end
        push!(bd, wall)
    end
    return Billiard(bd)
end

"""
    billiard_hexagonal_sinai(r, R, center = [0,0]; setting = "standard")
Create a sinai-like billiard, which is a hexagon of outer radius `R`, containing
at its center (given by `center`) a disk of radius `r`. The `setting` keyword
is passed to `billiard_polygon`.
"""
function billiard_hexagonal_sinai(r::Real = 0.5, R::Real = 1.0, center = [0,0];
    setting = "standard")
    r, R = promote(r, R)
    T = typeof(r); center = T[center...]
    bdr = billiard_polygon(6, R, center; setting = setting)
    DT = setting == "random" ? RandomDisk : Disk
    return Billiard(Disk(center, r), bdr...)
end




"""
    billiard_raysplitting_showcase(x=2.0, y=1.0, r1=0.3, r2=0.2) -> bd, rayspl
Showcase example billiard for ray-splitting processes. A rectangle `(x,y)` with a
SplitterWall at `x/2` and two disks at each side, with respective radii `r1`, `r2`.

**Notice**: This function returns a billiard `bd` as well as a `rayspl`
dictionary!
"""
function billiard_raysplitting_showcase(x=2.0, y=1.0, r1=0.3, r2=0.2)

    r1≥x/4 || r2≥x/4 && throw(ArgumentError("Disks overlap with walls!"))

    bdr =  billiard_rectangle(x, y)
    sw = SplitterWall([x/2, 0.0], [x/2,y], [-1,0], true)
    a1 = Antidot([x/4, y/2], r1, "Left Antidot")
    a2 = Antidot([3x/4, y/2], r2, "Right Antidot")

    sa = (φ, pflag, ω) -> pflag ? 0.5φ : 2.0φ
    Tp = (p) -> (φ, pflag, ω) -> begin
        if pflag
            p*exp(-(φ)^2/2(π/8)^2)
        else
            abs(φ) < π/4 ? (1-p)*exp(-(φ)^2/2(π/4)^2) : 0.0
        end
    end
    newoantidot = ((x, bool) -> bool ? -0.5x : -2.0x)
    newowall = ((x, bool) -> bool ? 0.5x : 2.0x)

    # I Need 3 ray-splitters because I want different affect for
    # different antidots
    raywall = RaySplitter([3], Tp(0.5), sa, newowall)
    raya = RaySplitter([1, 2], Tp(0.64), sa, newoantidot)

    return Billiard(a1, a2, sw, bdr...), (raywall, raya)
end


"""
    billiard_mushroom(l = 1.0, w = 0.2, r = 1.0, sloc = 0.0; door = true)
Create a mushroom billiard with cap radius `r`, stem width `w` and step
height `l`. The center of the cap (which is Semicircle) is always
at `[0, l]`. The center of the stem is located at `sloc`.

Optionally, the bottom-most `Wall` is a `Door` (see [`escapetime`](@ref)).
"""
function billiard_mushroom(stem_length = 1.0, stem_width=0.2, cap_radious=1.0,
    stem_location = 0.0; door = true)

    stloc = stem_location
    sl = stem_length; sw = stem_width; cr = cap_radious

    abs(stloc) + sw/2 > cr && error("Stem is outside the mushroom cap!")

    leftcorn = SV(-sw/2 + stloc, 0)
    rightcorn = SV(sw/2 + stloc, 0)
    upleftcorn = SV(-sw/2 + stloc, sl)
    uprightcorn = SV(sw/2 + stloc, sl)

    stembot = FiniteWall(leftcorn, rightcorn, SV(0, sw), door, "Stem bottom")
    stemleft = FiniteWall(upleftcorn, leftcorn, SV(sw, 0), false, "Stem left")
    stemright = FiniteWall(rightcorn, uprightcorn, SV(-sw, 0), false, "Stem right")


    farleft = SV(-cr, sl)
    farright = SV(cr, sl)

    capbotleft = FiniteWall(
    farleft, upleftcorn, SV(0, sw), false, "Cap bottom left")
    capbotright = FiniteWall(
    uprightcorn, farright, SV(0, sw), false, "Cap bottom right")

    cap = Semicircle([0.0, sl], cap_radious, [0.0, -1.0], "Mushroom cap")

    return Billiard(stembot, stemright, capbotright, cap, capbotleft, stemleft)
end

"""
    billiard_bunimovich(l=1.0, w=1.0)
Return a vector of `Obstacle`s that define a Buminovich billiard, also called a
stadium. The length is considered *without* the attached semicircles, meaning that the
full length of the billiard is `l + w`. The left and right edges of the stadium
are [`Semicircle`](@ref)s.

`billiard_stadium` is an alias of `billiard_bunimovich`.
"""
function billiard_bunimovich(l=1.0, w=1.0)

    l = convert(AbstractFloat, l)
    l, w = promote(l,w)
    o = typeof(l)(0.0)
    bw = InfiniteWall([o,o], [l,o], [o,  w], "Bottom wall")
    tw = InfiniteWall([l,w], [o,w], [o, -w], "Top wall")
    leftc = Semicircle([o, w/2], w/2, [l, o], "Left semicircle")
    rightc = Semicircle([l, w/2], w/2, [-l, o], "Right semicircle")
    return Billiard(bw, rightc, tw, leftc)
end

billiard_stadium = billiard_bunimovich
