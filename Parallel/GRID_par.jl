using Base.Iterators, Gtk, LinearAlgebra, Plots, BenchmarkTools, MatrixMarket
pyplot()

function rmin(a,b)
    if a < b
        return a
    else 
        return b
    end
end

function rmax(a,b)
    if a < b
        return b
    else 
        return a
    end
end

"""
Generating a matrix of dimension n.
"""

A = [3-1im 1+5im 4 1 ; 0 -2+3im -7im 8-2im; 1 10-5im 0 5; 2 8 2im -1im]

function generate_matrix(a, b, n)
    Ar =  a .+ (b-a) * rand(n,n)
    A = fill(1.0im, n, n)
    for i in 1:n 
        for j in 1:n
            A[i,j] = Ar[i,j] + (a + (b-a) * rand()) * 1.0im
        end
    end
    return A
end
"""
Building a grid of grid_size^2 elements that will encompass the pseudospectra of perturbance eps of matrix A
according to Gershgorin's theorem.
"""

function generate_grid(A, eps)
    max_cord = 0
    min_cord = 0
    temp_radius = 0
    radius = 0


   #Calculating the radius R'
    for i in 1:size(A, 1)
        for j in 1:size(A, 2)
            if i!=j
                temp_radius = 0
                temp_radius += abs(A[i,j])
                radius = max(radius, temp_radius)
            end
        end
    end
    radius+= sqrt(size(A,1))*eps

    #Determining min and max limits value of the square grid
    DiagA = LinearAlgebra.Diagonal(A)
    egv = eigvals(DiagA)

    min_real_cord = minimum(real(egv)) - radius
    min_imag_cord = minimum(imag(egv)) - radius
    
    min_cord = floor(rmin(min_real_cord, min_imag_cord))


    max_real_cord = maximum(real(egv)) + radius
    max_imag_cord = maximum(imag(egv)) + radius

    max_cord = floor(rmax(max_real_cord, max_imag_cord))+1


    grid_size = 100
    X = range(min_cord, max_cord, grid_size)
    Y = range(max_cord, min_cord, grid_size)

    return X, Y
end

"""
Calculate sigmin
"""
function sigmin_cpt(A)
    svd = svdvals(A)
    return svd[length(svd)]
end

"""
Calculate sigmax
"""

function sigmax_cpt(A)
    svd = svdvals(A)
    return svd[1]
end

"""
GRID Algorithm : Calculating minimal singular value for every point of the grid.
"""
function svd_compute(A, eps)
    (X, Y) = generate_grid(A, eps)
	sigmin = zeros(length(X), length(Y))
    Threads.@threads for i in 1:100
        for j in 1:100
            sigmin[j,i] = sigmin_cpt((X[i] + Y[j] * 1im)*I - A)
        end
    end
	return sigmin, X, Y
end


"""
plotting the results of GRID algo using contour and scatter.
"""


function contours(A, eps)
    (z, X, Y) = svd_compute(A, eps)
	data = Plots.contour(X, Y, z, colors = "grey", levels = [0, eps], fill = false, colorbar = false)
    scatter!(eigvals(A), lab="eigenvalues")
	plot = plot!(data, lab="pseudospectra", fmt= :png)
    title!(string("The Pseudospectra of Matrix A, ε = ", string(eps)))
    xaxis!("Re")
    yaxis!("Im")
    filename = string("contour", string(eps),".png")
    Plots.savefig(plot, filename)
    return filename
end


"""
Showing the results in a GUI.
"""
function execute(A)
    pswindow = GtkWindow("GRID algorithm", 400, 100)
    grid = GtkGrid()
    
    label = GtkLabel("label")
    GAccessor.markup(label, "<b>Please enter the perturbance ε value in the textbox.</b>")  

    eps_entry = GtkEntry()

    grid[1,2] = label
    grid[1,4] = eps_entry

    set_gtk_property!(grid, :column_homogeneous, true)
    set_gtk_property!(grid, :column_spacing, 15)

    push!(pswindow, grid)
    showall(pswindow)
    
    signal_connect(eps_entry, "key-press-event") do widget, event
        if event.keyval == 65293 #Enter
            plot_window = GtkWindow("Resulting", 400, 300)
            eps = 0
            eps = parse(Float64, get_gtk_property(eps_entry,:text,String))
            set_gtk_property!(eps_entry, :text, "")
            image = GtkImage(contours(A, eps))
            push!(plot_window, image)
            show(image)
            showall(plot_window)
        end
    end
end



"""
TESTS

A = MatrixMarket.mmread("qc324.mtx")
A = Matrix(A)
A = A[1:10, 1:10]
@btime svd_compute(A, 1)"""