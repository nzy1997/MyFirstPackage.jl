"""
    D3Q19 <: AbstractLBConfig{3, 19}

A lattice Boltzmann configuration for 3D, 19-velocity model.
"""
struct D3Q19 <: AbstractLBConfig{3, 19} end

directions(::D3Q19) = (
        Point(1, 0, -1),
        Point(0, 1, 1),
        Point(0, 1, -1),
        Point(1, 0, 1), 
        Point(1, -1, 0),
        Point(1, 1, 0),
        Point(0, 0, 1),
        Point(0, 1, 0),
        Point(1, 0, 0),
        Point(0, 0, 0), 
        Point(-1, 0, 0),
        Point(0, -1, 0),
        Point(0, 0, -1),
        Point(-1, -1, 0),
        Point(-1, 1, 0),
        Point(-1, 0, -1),
        Point(0, -1, 1),
        Point(0, -1, -1),
        Point(-1, 0, 1)
        )

# directions[k] is the opposite of directions[flip_direction_index(k)
function flip_direction_index(::D3Q19, i::Int)
    return 20 - i
end

weights(::D3Q19) = (1/52,1/52,1/52,1/52,1/52,1/52,1/13,1/13,1/13,4/13,1/13,1/13,1/13,1/52,1/52,1/52,1/52,1/52,1/52)

"""
	LatticeBoltzmann{D, N, T, CFG, MT, BT}

A lattice Boltzmann simulation with D dimensions, N velocities, and lattice configuration CFG.
"""
struct LatticeBoltzmann3d{D, N, T, CFG <: AbstractLBConfig{D, N}, MT <: AbstractArray{Cell{N, T}}, BT <: AbstractArray{Bool}}
	config::CFG # lattice configuration
	grid::MT    # density of the fluid
	gridcache::MT # cache for the density of the fluid
	barrier::BT # barrier configuration
end

function LatticeBoltzmann3d(config::AbstractLBConfig{D, N}, grid::AbstractArray{<:Cell}, barrier::AbstractArray{Bool}) where {D, N}
	@assert size(grid) == size(barrier)
	return LatticeBoltzmann3d(config, grid, similar(grid), barrier)
end


function stream!(
    lb::AbstractLBConfig{3, N},  # lattice configuration for 3D
    newgrid::AbstractArray{D}, # the updated grid (3D)
    grid::AbstractArray{D}, # the original grid (3D)
    barrier::AbstractArray{Bool} # the barrier configuration (3D)
) where {N, T, D<:Cell{N, T}}
ds = directions(lb)  # Ensure this provides the D3Q19 velocity vectors
@inbounds for ci in CartesianIndices(newgrid)
    i, j, k = ci.I  # Now we have three indices for 3D
    newgrid[ci] = Cell(ntuple(N) do d  # Collect the densities
        ei = ds[d]
        l, m, n = size(grid)  # Dimensions for 3D grid
        i2, j2, k2 = mod1(i - ei[1], l), mod1(j - ei[2], m), mod1(k - ei[3], n)
        if barrier[i2, j2, k2]
            # if the cell is a barrier, the fluid flows back
            density(grid[i, j, k], flip_direction_index(lb, d))
        else
            # otherwise, the fluid flows to the neighboring cell
            density(grid[i2, j2, k2], d)
        end
    end)
end
end

function step!(lb::LatticeBoltzmann3d)
	copyto!(lb.gridcache, lb.grid)
	stream!(lb.config, lb.grid, lb.gridcache, lb.barrier)
	lb.grid .= collide.(Ref(lb.config), lb.grid)
	return lb
end
