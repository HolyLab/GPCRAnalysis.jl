function markers(modelnum::Integer, positions::AbstractVector{<:AbstractVector{<:Real}}, radius::Real, color)
    Base.depwarn("`markers` is deprecated, use `marker.(modelnum, positions, radius, color)` instead", :markers)
    return [marker(modelnum, pos, radius, color) for pos in positions]
end
