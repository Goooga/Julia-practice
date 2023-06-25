include("practice1.jl")
using LinearAlgebra

#1
function pow(a::T, n::M) where {T, M}
    res = one(T)
    while n > 0
        if mod(n, 2) == 1
            res *= a
        end
        a*=a
        n = div(n, 2)
    end
    return res
end

#2
Base. one(::Type{Matrix{T}}) where{T} = [1 0; 0 1]

function fibonacci(n)
    m = [ 1 1; 1 0]
    return pow(m, n - 1)[1, 1]
end

#3
function log_(a::T,x::T,eps=0.001) where T<:Union{AbstractFloat, Rational}
    if a==1 || a < 0
        return nothing
    end
    if a<1
        return -log_(1/a,x)
    end
    @assert x > 0

    z=x
    t=one(T)
    y=zero(T)
    #ИНВАРИАНТ:  x = z^t * a^y
    while z < 1/a || z > a || t > eps 
        if z<1/a
            z*=a
            y-=t
        elseif z>a
            z/=a
            y+=t
        elseif t>eps
            t/=2
            z*=z
        end
    end
    return y 
end

#4
function bisection(f::Function, a, b, epsilon)
    @assert f(a)*f(b) < 0 
    @assert a < b
    f_a = f(a)
    #ИНВАРИАНТ: f_a*f(b) < 0
    while b-a > epsilon
        t = (a+b)/2
        f_t = f(t)
        if f_t == 0
            return t
        elseif f_a*f_t < 0
            b=t
        else
            a, f_a = t, f_t
        end
    end
    return (a+b)/2
end

#5
# println(bisection(x->cos(x),-pi,pi/4,0.0001))

#6
function newton(r::Function, x, epsilon, num_max = 50)
    k=0
    dx = -r(x)
    while abs(dx) > epsilon && k <= num_max
        dx = -r(x)
        x += dx
        k += 1
    end
    k > num_max && @warn("Требуемая точность не достигнута")
    return x
end

#7
# println(newton(x->cos(x), pi/4, 0.001 ))

#8
# p1 = make_polynom([2, 3, 5])
# println(display(p1))
# println(newton(x->/(polyval(p1, x)...),1,0.0001))
# println(newton(x->p1(x),1,0.001))

Base.@kwdef struct NewtonFractal
    basin_colors::NTuple{N, Symbol} where N  
        # - при создании объекта требуется задать цвета раскраски бассейнов притяжений, 
        # например, NewtonFractal((:red, :green, :blue))
    order = length(basin_colors) # порядок уравнения
    roots = tuple((exp(im*2π/order*i) for i ∈ 0:order-1)...) # корни уравнения: z^order = 1
    square_side = 4 # размер квадрата
    points_number = 1_000_000
    ε = 0.3 # определяет окрестность "притяжения" корня
    iterations_number = 10 
        # - при построении фрактала после этого числа итераций определяется, в окрестность 
        # неотвратимого притяжения какого корня попала точка (или не попала ни в какую окрестность)
    basin = tuple((Complex{Float64}[] for _ in 1:order)...)
end 

function build!(fractal::NewtonFractal)
    @assert fractal.order > 1
    n = fractal.order
    r(z) = (z-1/z^(n-1))/n # r = f(z)/f'(z), where f(z) = z^n-1

    for _ in 1:fractal.points_number
        z_0 = fractal.square_side*(rand(Complex{Float64}) - Complex{Float64}(0.5, 0.5))
        z = z_0
        for _ in 1:fractal.iterations_number
            z -= r(z)
        end    
        i = findfirst(abs.(fractal.roots .- z) .<= fractal.ε)
        !isnothing(i) && push!(fractal.basin[i], z_0)
    end
end

# using GLMakie
# function draw(fractal::NewtonFractal)
#     f = Figure()
#     a = Axis(f[1,1])
#     for i ∈ 1:fractal.order
#         scatter!(a, reim(fractal.basin[i])..., color = fractal.basin_colors[i], markersize = 1)
#     end
#     a.aspect = 1
#     display(f)
# end

# #---------------ТЕСТ-----------------------
# fractal = NewtonFractal(
#     basin_colors=(:red,:green,:blue),
#     iterations_number = 20, 
#     ε=0.5,
#     points_number = 10_000_000
# )
# build!(fractal)
# draw(fractal)


# r(x) = x-> /(polyval(p,x)...)