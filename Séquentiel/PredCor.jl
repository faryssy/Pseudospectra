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
    print(Ar)
    A = fill(1.0im, n, n)
    for i in 1:n 
        for j in 1:n
            A[i,j] = Ar[i,j] + (a + (b-a) * rand()) * 1.0im
        end
    end
    return A
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
Prediction-Correction Method for a pseudospectra section
"""


function pred_cor(A, eps, l0, k, d0, tol)

    
    Z  = fill(0.0 + 1.0im, k)
    t0 = eps
    z_new = l0 + t0 * d0
    s_min = 0
    v_min = fill(0 + 1*im, size(A, 2))
    u_min = fill(0 + 1*im, size(A, 2))

    while abs(sigmin_cpt(z_new*I - A) - eps) > tol * eps
        z_old = z_new
        M = z_old*I - A
        (U, S, V) = svd(M)
        s_min = S[length(S)]
        u_min = U[:, size(M, 1)]
        v_min = V[:, size(M, 1)]
        d0 = v_min' * u_min
        z_new = z_old - (s_min - eps)/real(conj(d0)*v_min'*u_min)*d0

    end

    Z[1] = z_new

    for i in 2:k
        # Prediction
        r = (1im * v_min'* u_min)/abs(v_min' * u_min)
        tau = 0.1
        zt = Z[i-1] + tau * r

        # Correction
        M = zt*I - A
        (U, S, V) = svd(M)
        s_min = S[length(S)]
        u_min = U[:, size(M, 1)]
        v_min = V[:, size(M, 1)]
        Z[i] = zt - (s_min - eps)/(u_min' * v_min)    
    end

    return Z
end

function pred_cor_global(A, eps, k)
    egv = eigvals(A)
    Ps = fill(0.0 + 1.0im, length(egv)*k) 
    Z = fill(0.0 + 1.0im, k)
    for i in 0 : length(egv)-1
        Z = pred_cor(A, eps, egv[i+1], k, 1, 0.001)
        Ps[i*k+1:(i+1)*k] = Z
    end
    return Ps
end

"""
Plotting the results of Prediction-Correction algo
"""
function plots(A, eps, k)
    egv = eigvals(A)
    plotting = nothing
    plot()
    Z = pred_cor_global(A, eps, k)
    for i in 0:length(egv)-1
        Re_Z = real(Z[i*k+1:(i+1)*k])
        Im_Z = imag(Z[i*k+1:(i+1)*k])
        data = Plots.plot!(Re_Z, Im_Z, label = "pseudospectra", color ="black")
        plotting = plot!(data, fmt= :png) 
    end
    scatter!(eigvals(A), lab="eigenvalues", color="red")
    title!(string("Pred-Cor : Pseudospectra of Matrix A, perturbance = ", string(eps)))
    xaxis!("Re")
    yaxis!("Im")
    filename = string("plot", string(eps),".png")
    Plots.savefig(plotting, filename)
    return filename
end
"""
Showing the results in a GUI.
"""
function execute(A, k)
    pswindow = GtkWindow("Pred-cor algorithm", 400, 100)
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
            image = GtkImage(plots(A, eps, k))
            push!(plot_window, image)
            show(image)
            showall(plot_window)
        end
    end
end

    

plots(A, 1, 295) 
#execute(A)
#A = generate_matrix(-10, 10, 100)

"""
TESTS
"""
A = MatrixMarket.mmread("qc324.mtx")
A = Matrix(A)
A = A[1:20, 1:20]
@btime pred_cor_global(A, 1, 25)