"""
    project_to_line(point, c, n)
Project given `point` to line that contains point `c` and has **normal vector** `n`.
"""
function project_to_line(point, c, n)
    posdot = dot(c - point, n)
    intersection = point + posdot*n
end
