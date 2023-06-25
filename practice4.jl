include("practice2.jl")
# using GLMakie
using LinearAlgebra
using Plots
# 1
function exp_macloren_n_accuracy(x, n)
    res = 1
    term = 1
    for i in 1:n
        term = term * x / i
        res += term
    end
    return res
end

# println(exp_macloren_n_accuracy(2, 2))
#2
function exp_(x)
    if x == trunc(x)
        return pow(ℯ, x)
    end
    res = 1
    term = 1
    i = 1
    while term > eps()
        term = term*x/i
        res += term
        i+=1
    end
    return res
end
println(exp_(2))

# 3
function bessel(x, order)
    i = 1
    sum = 0
    term = 1/factorial(order)
    while sum + term != sum
        sum += term
        term = -term * (x*x) /(i * (order + i)*4)
        i+=1
    end
    return sum *(x/2)^order
end

values = 0:0.1:20
p = plot()
for order in 0:5
	plot!(p, values, bessel.(values, order))
end
display(p)
# 4
function reverse_gauss_jordan(matrix::AbstractMatrix{T}, b::AbstractVector{T} ) where T
    x = similar(b)
    n = size(matrix, 1) #количество строк
	for i in 0:n-1
		x[n-i] = (b[n-i] - sum(matrix[n-i,n-i+1:end] .* x[n-i+1:end])) / matrix[n-i,n-i]
	end
    return x
end

function swap_columns(matrix::AbstractMatrix{T}, index1, index2) where T
    rows_count = size(matrix, 1)
    for i in 1:rows_count
        matrix[i, index1], matrix[i, index2] = matrix[i, index2], matrix[i, index1]
    end
end 

# 5
function stepped_view_matrix(matrix::AbstractMatrix{T}, b::AbstractVector{T} = zeros(T, size(matrix, 1))) where T
    matrix2 = transpose(matrix)
    b = transpose(b)
   for j in 1:size(matrix2, 2) - 1
        if matrix2[j, j] == 0
            for k in j+1:size(matrix2, 2)
                columns_is_swaped = false
                if matrix2[j, k] !=  0
                    swap_columns(matrix2, k, j)
                    swap_columns(b, k, j)
                    columns_is_swaped = true
                end
            end
            if columns_is_swaped == false
                break
            end
        end
        for i in j+1:size(matrix2, 2)
            matrix2[j:end, i] -= @view(matrix2[j:end, j]) * matrix2[j, i]/matrix2[j,j]
            b[i] -= b[j] * matrix2[j, i]/ matrix2[j, j]
        end
    end
    return transpose(matrix2), transpose(b)
end

@inline function sumprod(vec1::AbstractVector{T}, vec2::AbstractVector{T})::T where T
	s = zero(T)
	@inbounds for i in eachindex(vec1)
	s = fma(vec1[i], vec2[i], s) # fma(x, y, z) вычисляет выражение x*y+z
	end
	return s
end

#6
function jordan_gauss(matrix::AbstractMatrix{T}, b::AbstractVector{T} = zeros(T,size(matrix, 1))) where T
    #сначала приводим матрицу к ступенчатому виду
    tup = stepped_view_matrix(matrix, b)
    
    x = tup[2]
    n = size(matrix, 1)
    for i in 0:n-1
		x[n-i] = (b[n-i] - sumprod(@view(matrix[n-i,n-i+1:end]), @view(x[n-i+1:end]))) / matrix[n-i,n-i]
	end

	return x
end

#7.
# for n in 50:50:1000
# 	println("Матрица порядка ",n,"×",n,":")
# 	@time jordan_gauss(randn(n,n),randn(n))
# 	Матрица порядка 50×50:
#   0.001229 seconds (4.90 k allocations: 1.607 MiB)
# Матрица порядка 100×100:
#   0.024207 seconds (19.80 k allocations: 11.768 MiB, 63.37% gc time)
# Матрица порядка 150×150:
#   0.037570 seconds (44.70 k allocations: 38.288 MiB, 27.63% gc time)
# Матрица порядка 200×200:
#   0.079702 seconds (79.60 k allocations: 89.536 MiB, 20.95% gc time)
# Матрица порядка 250×250:
#   0.100605 seconds (124.50 k allocations: 173.552 MiB, 17.84% gc time)
# Матрица порядка 300×300:
#   0.198107 seconds (179.40 k allocations: 297.661 MiB, 15.62% gc time)
# Матрица порядка 350×350:
#   0.235374 seconds (244.30 k allocations: 469.036 MiB, 18.55% gc time)
# Матрица порядка 400×400:
#   0.329269 seconds (319.20 k allocations: 695.288 MiB, 18.09% gc time)
# Матрица порядка 450×450:
#   0.416443 seconds (404.10 k allocations: 984.247 MiB, 17.57% gc time)
# Матрица порядка 500×500:
#   0.905616 seconds (499.00 k allocations: 1.312 GiB, 29.62% gc time)
# Матрица порядка 550×550:
#   0.797382 seconds (603.90 k allocations: 1.738 GiB, 16.46% gc time)
# Матрица порядка 600×600:
#   0.971537 seconds (718.80 k allocations: 2.248 GiB, 14.88% gc time)
# Матрица порядка 650×650:
#   1.333013 seconds (843.70 k allocations: 2.850 GiB, 24.59% gc time)
# Матрица порядка 700×700:
#   1.777913 seconds (978.60 k allocations: 3.549 GiB, 12.58% gc time)
# Матрица порядка 750×750:
#   2.467444 seconds (1.12 M allocations: 4.354 GiB, 11.68% gc time)
# Матрица порядка 800×800:
#   2.464682 seconds (1.28 M allocations: 5.272 GiB, 19.32% gc time)
# Матрица порядка 850×850:
#   2.549380 seconds (1.44 M allocations: 6.312 GiB, 11.87% gc time)
# Матрица порядка 900×900:
#   3.220300 seconds (1.62 M allocations: 7.479 GiB, 15.54% gc time)
# Матрица порядка 950×950:
#   3.391912 seconds (1.80 M allocations: 8.782 GiB, 11.44% gc time)
# Матрица порядка 1000×1000:
#   3.969842 seconds (2.00 M allocations: 10.227 GiB, 14.31% gc time)
# end

# 8 
function rang(matrix::AbstractMatrix{T}) where T
    tup = stepped_view_matrix(copy(matrix))
    matrix2 = tup[1]
    k = 0
    for j in 1:size(matrix2, 2)
        iszeros = true
        for i in matrix2[1:end, j]
            if i != 0
                iszeros = false
                break
            end
        end
        if iszeros == false
            k+=1
        end
    end
    return k
end

#9
function det(matrix::AbstractMatrix{T}) where T
    stepped_view_matrix(matrix)

    det = oneunit(T)
    i = 1

    while i <= size(matrix, 1)
		if matrix[i, i] == zero(T)
			return 0
		end
		det *= matrix[i, i]
		
		i += 1
    end

    return det
end

a = [3. 4. 5. 6; 7. 8. 9. 10; 11. 12. 13. 14.; 14. 15. 16. 0.]
println(stepped_view_matrix(a)[1])
println(jordan_gauss(a, [1., 1., 1., 1.]))
println(rang(a))
m = [-2. 3.; 1. 4.]
println(det(m))