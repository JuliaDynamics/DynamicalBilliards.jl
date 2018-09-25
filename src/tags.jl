module Tagging
export tag
using DynamicalBilliards
CONCAVE = Union{Semicircle}

tag(p::Particle) = "[STRAIGHT]"
tag(p::MagneticParticle) = "[MAGNETIC]"

function tag(bd::Billiard)
    s = DynamicalBilliards.isperiodic(bd) ? "[PERIODIC]" : ""
    if length(bd) > 8
        s *= "[OMNI]"
    elseif any(o -> typeof(o) <: CONCAVE, bd)
        s *= "[CONCAVE]"
    end
    return s
end

tag(t::Union{<: RaySplitter, <: Tuple}) = "[RAYSPLIT]"

end
