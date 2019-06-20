module HoloJulia

using Images

export HoloCrea, HoloMemberPC, HoloEnsemblePC

"""
     HoloCrea(U, lines, alf, kind)
Creates the hologram of the desired field U """
function HoloCrea(U, lines, alf; kind=1)

    # SPATIAL SCALE
    Ux=15.36E-3;
    Uy=8.64E-3;
    Hsize=size(U,2);
    Vsize=size(U,1);

    # Generates ranges for xs and ys
    ys = Ux*(2/Hsize)*collect(range(-Hsize/2,length=Hsize,stop=Hsize/2-1))
    xs = Uy*(2/Vsize)*collect(range(-Vsize/2,length=Vsize,stop=Vsize/2-1))

    # Coordinates
    mx, ny = length(xs), length(ys)
    Ys = reshape(xs, mx, 1)
    Xs = reshape(ys, 1, ny)

    # Hologram creation: receives complex field U properly normalized such that U is in [0,1]
    A = abs.(U)
    g = A./maximum(A)
    phase = angle.(U)

    # lines = 100
    linespMM = 2*lines/(1E-2)
    # alf = 0
    kxs = linespMM*cos(alf)
    kys = linespMM*sin(alf)
    if kind==0
        H = 1 .*(mod.(phase .+ (kxs*Xs .+ kys*Ys), 2*pi)./pi .- 1) .+ 1;
    elseif kind==1
        H = g .*(mod.(phase .+ (kxs*Xs .+ kys*Ys), 2*pi)./pi .- 1) .+ 1;
    elseif kind==2
        H = 1 .*(mod.(phase .+ (kxs*Xs .+ kys*Ys), 2*pi)./pi .- 1) .+ 1;
        H = round.(H./maximum(H));
    else
        H = g .*(mod.(phase .+ (kxs*Xs .+ kys*Ys), 2*pi));
        H = round.(H./maximum(H));
    end
    H = H./maximum(H);

    return H
end

"""
     HolEnsemblePC(Useed, Vsize, Hsize, XYsize, circle, N, Ne, lines, alf, dirname, filename, kk)
Creates many holograms of the desired randomly displaced field U """
function HoloEnsemblePC(Useed, Vsize, Hsize, XYsize, circle, N, Ne, lines, alf, dirname, filename, kk)
    # Generates the desired holograms according to the given parameters
    for memberE = 1:Ne
        # Computes the field of one member
        Upc = HoloMemberPC(Useed, Vsize, Hsize, XYsize, circle, N)

        # Hologram creation
        Hologram = HoloCrea(Upc, lines, alf, kind=kk)
        save("$dirname/$(filename)-$(lpad(memberE,3,'0')).jpg", colorview(Gray, Hologram))
    end
end

"""
     memberPC(Useed, Xs, Ys, w0, TopologicalCharge, RadialOrder, circle, N)

Simulates a member of the ensamble to generate a partially coherent Useed beam """
function HoloMemberPC(Useed, ren, col, XYsize, circle::Float64, N::Int64)

# Preallocates memory for matrices U and Uaux
Uaux = zeros(ComplexF64, ren, col)
U = zeros(ComplexF64, ren, col)

# Uniformly distributed vortices in a circle of radius "circle"
r = circle*sqrt.(rand(N,1))
theta = 2*pi*rand(N,1)
phi = 2*pi*rand(N,1)
xv = Integer.(div.(r.*cos.(theta),XYsize))
yv = Integer.(div.(r.*sin.(theta),XYsize))

# Coherent superposition
for member=1:N
    Uaux = exp(im * phi[member]) .* MatrixTranslate(Useed, xv[member], yv[member])
    Uaux = Uaux./maximum(abs.(Uaux))
    U = Uaux + U
end

return U
end

"""
    MatrixTranslate(U, xv::Int64, yv::Int64)

Dirty matrix translation """
function MatrixTranslate(U, xv::Int64, yv::Int64)

# Preallocates memory for matrices U and Uaux
ren = size(U,1)
col = size(U,2)
V = zeros(ComplexF64, ren, col)

# Shifts!
if xv>=0 && yv>=0
    V[yv+1:end, xv+1:end] = U[1:end-yv, 1:end-xv]

elseif xv<=0 && yv>=0
    V[yv+1:end, 1:end-abs(xv)] = U[1:end-yv, abs(xv)+1:end]

elseif xv>=0 && yv<=0
    V[1:end-abs(yv), xv+1:end] = U[abs(yv)+1:end, 1:end-xv]

elseif xv<=0 && yv<=0
    V[1:end-abs(yv), 1:end-abs(xv)] = U[abs(yv)+1:end, abs(xv)+1:end]
end

return V

end


end # module
