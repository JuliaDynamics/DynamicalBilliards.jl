obcolor(::Obstacle) = (0,0.6,0)
obcolor(::Union{RandomWall, RandomDisk}) = (149/255, 88/255, 178/255)
obcolor(::Union{SplitterWall, Antidot, Ellipse}) = (0.8,0.0,0)
obcolor(::PeriodicWall) = (0.8,0.8,0)
obalpha(::Obstacle) = 0.5
obalpha(::Union{Antidot, Ellipse}) = 0.1
obls(::Obstacle) = "solid"
obls(::Union{SplitterWall, Antidot, Ellipse}) = "dashed"
obls(::PeriodicWall) = "dotted"

"""
    plot(obst::Obstacle; kwargs...)
Plot given obstacle on the current `PyPlot` axes.

The default arguments for each type of obstacle have been chosen for maximum
clarity and consistency.

The `kwargs...` given by the user are keywords passed directly into PyPlot's
constructors. For `Wall` obstacles, kwargs are passed into `PyPlot.plot()`. For
`Circular` obstacles, kwargs are passed into `matplotlib.patches.Circle` or `Arc`.
"""
function plot(d::Obstacle) end

function plot(d::Circular; kwargs...)
    edgecolor = obcolor(d)
    facecolor = (edgecolor..., obalpha(d))
    circle1 = PyPlot.plt.Circle(d.c, d.r;
        edgecolor = edgecolor, facecolor = facecolor,
        linestyle = obls(d), lw = 2.0, kwargs...)
    PyPlot.gca().add_artist(circle1)
end

function plot(d::Semicircle; kwargs...)
    theta1 = atan(d.facedir[2], d.facedir[1])*180/π + 90
    theta2 = theta1 + 180
    edgecolor = obcolor(d)
    s1 = PyPlot.matplotlib.patches.Arc(d.c, 2d.r, 2d.r, theta1 = theta1, theta2 = theta2, edgecolor = edgecolor,
    lw = 2.0, kwargs...)
    PyPlot.gca().add_artist(s1)
end

function plot(w::Wall; kwargs...)
    if typeof(w) <: FiniteWall &&  w.isdoor
       PyPlot.plot([w.sp[1],w.ep[1]],[w.sp[2],w.ep[2]];
       color="black", linestyle = "-", lw = 2.0, kwargs...)
       PyPlot.plot([w.sp[1],w.ep[1]],[w.sp[2],w.ep[2]];
       color=(0, 0.9, 0.9), linestyle = "--", lw = 2.0, kwargs...)
    else
        PyPlot.plot([w.sp[1],w.ep[1]],[w.sp[2],w.ep[2]];
        color=obcolor(w),
        linestyle = obls(w), lw = 2.0, kwargs...)
    end
end

function plot(e::Ellipse; kwargs...)
    edgecolor = obcolor(e)
    facecolor = (edgecolor..., obalpha(e))
    ellipse = PyPlot.matplotlib.patches.Ellipse(e.c, 2e.a, 2e.b;
        edgecolor = edgecolor, facecolor = facecolor,
        linestyle = obls(e), lw = 2.0, kwargs...)
    PyPlot.gca().add_artist(ellipse)
end
