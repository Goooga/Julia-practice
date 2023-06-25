# include("Plots.jl")

using Plots

using LinearAlgebra

# 1.

Vector2D{T <: Real} = NamedTuple{(:x, :y), Tuple{T,T}}

Base. +(a::Vector2D{T},b::Vector2D{T}) where T = Vector2D{T}(Tuple(a) .+ Tuple(b))
Base. -(a::Vector2D{T}, b::Vector2D{T}) where T = Vector2D{T}(Tuple(a) .- Tuple(b))
Base. *(α::T, a::Vector2D{T}) where T = Vector2D{T}(α.*Tuple(a))

# norm(a) - длина вектора, эта функция опредедена в LinearAlgebra
LinearAlgebra.norm(a::Vector2D) = norm(Tuple(a))

# dot(a,b)=|a||b|cos(a,b) - скалярное произведение, эта функция определена в LinearAlgebra
LinearAlgebra.dot(a::Vector2D{T}, b::Vector2D{T}) where T = dot(Tuple(a), Tuple(b))

Base. cos(a::Vector2D{T}, b::Vector2D{T}) where T = dot(a,b)/norm(a)/norm(b)

# xdot(a,b)=|a||b|sin(a,b) - косое произведение
xdot(a::Vector2D{T}, b::Vector2D{T}) where T = a.x*b.y-a.y*b.x

Base.sin(a::Vector2D{T}, b::Vector2D{T}) where T = xdot(a,b)/norm(a)/norm(b)
Base.angle(a::Vector2D{T}, b::Vector2D{T}) where T = atan(sin(a,b),cos(a,b))
Base.sign(a::Vector2D{T}, b::Vector2D{T}) where T = sign(sin(a,b))

Segment2D{T <: Real} = NamedTuple{(:A, :B), NTuple{2,Vector2D{T}}}

# 2
function oneside(P::Vector2D{T}, Q::Vector2D{T}, s::Segment2D{T})::Bool where T
	l = s.B - s.A #направление

	return sin(l, P-s.A) * sin(l,Q-s.A) > 0
end

# 3. 
oneside(F::Function, P::Vector2D, Q::Vector2D)::Bool = ( F(P...) * F(Q...) > 0 )

# 4.
isinner(P::Vector2D, s::Segment2D)::Bool = (s.A.x <= P.x <= s.B.x || s.A.x >= P.x >= s.B.x) && (s.A.y <= P.y <= s.B.y || s.A.y >= P.y >= s.B.y)

function intersection(s1::Segment2D{T},s2::Segment2D{T})::Union{Vector2D{T},Nothing} where T
	A = [s1.B[2]-s1.A[2] s1.A[1]-s1.B[1]
		s2.B[2]-s2.A[2] s2.A[1]-s2.B[1]]

	b = [s1.A[2]*(s1.A[1]-s1.B[1]) + s1.A[1]*(s1.B[2]-s1.A[2])
		s2.A[2]*(s2.A[1]-s2.B[1]) + s2.A[1]*(s2.B[2]-s2.A[2])]

	x,y = A\b

	# Замечание: Если матрица A - вырожденная, то произойдет ошибка времени выполнения
	if isinner((;x, y), s1)==false || isinner((;x, y), s2)==false
		return nothing
	end

	return (;x, y) #Vector2D{T}((x,y))
end

# println(intersection((A=(x=1.0, y=0.0), B=(x=0.0, y=1.0)), (A=(x=0.0, y=-1.), B=(x=1.0, y=1.0))))
# 5. Написать функцию, проверяющую лежит ли заданная точка внутри заданного многоугольника.

function is_dot_in_polygon(point::Vector2D{T}, polygon::AbstractArray{Vector2D{T}}) where T

	if length(polygon) < 3
        return false
    end

	angle_sum = zero(Float64)

	# Вычислить алгебраическую (т.е. с учетом знака) сумму углов, между направлениями из заданной точки на каждые две сосоедние вершины многоугольника.

	for i in firstindex(polygon):lastindex(polygon)
		angle_sum += angle(polygon[i] - point, polygon[i % lastindex(polygon) + 1] - point) 
	end
	
	return abs(angle_sum) > π # поскольку эта сумма не может иметь промежуточных значений между указанными величинами, то сравнивать можно с значением pi)
end

println(is_dot_in_polygon( (x=5,y=5), [(x = 0, y = 0), (x = 4, y = 0), (x = 4, y = 4), (x = 0, y = 4)]))

# 6. Написать функцию, проверяющую, является ли заданный многоугольник выпуклым.

function is_convex_polygon(coords::AbstractArray{Vector2D{T}}) where T

    n = length(coords)

    if n < 3
        return false
    end
    
    signs = zeros(n) # Считаем знаки всех углов многоугольника

    for i in 1:n
        p1 = coords[i]
        p2 = coords[mod1(i+1, n)]
        p3 = coords[mod1(i+2, n)]
        cross_prod = (p2.x - p1.x) * (p3.y - p2.y) - (p3.x - p2.x) * (p2.y - p1.y)
        signs[i] = sign(cross_prod)
    end
	for i in 2:n
		if signs[i] != signs[1]
			return false
		end
	end
    return true
end

coords = [
    (x = -1, y = -1),
    (x = 2, y = 0),
    (x = 1, y = 1),
    (x = 2, y = 2),
    (x = 0, y = 2)
]

println(is_convex_polygon(coords))

# 7. Выпуклая оболочка по Джарвису 
function next!(convex_shell2::AbstractVector{Int64}, points2::AbstractVector{Vector2D{T}}, ort_base::Vector2D{T})::Int64 where T
	cos_max = typemin(T)
	i_base = convex_shell2[end] # Получение индекса последней точки, которая будет использоваться в цикле для сравнения с другими точками.
	resize!(convex_shell2, length(convex_shell2) + 1)
	for i in eachindex(points2) # Проходимся по всем точкам
		if points2[i] == points2[i_base] # из points2[i] == points2[i_base] не следует что i == i_base
			continue
		end
		ort_i = points2[i] - points2[i_base] # - не нулевой вектор, задающий направление на очередную точку
		cos_i = cos(ort_base, ort_i) # косинус между вектором, определяющим направление от базовой точки к текущей точке
		if cos_i > cos_max           # и вектором, определяющим текущее направление для поиска следующей точки.
			cos_max = cos_i
			convex_shell2[end] = i
		elseif cos_i == cos_max && dot(ort_i, ort_i) > dot(ort_base, ort_base) # на луче, содержащем сторону выпуклого многоугольника, может оказаться более двух точек заданного множества (надо выбрать самую дальнюю из них)
			convex_shell2[end] = i
		end
	end
	return convex_shell2[end]
end

function jarvis!(dots::AbstractArray{Vector2D{T}})::AbstractArray{Vector2D{T}} where T

    @assert length(dots) > 1
	
    ydata = [dots[i].y for i in 1:length(dots)]
    start_dot = findmin(ydata) # стартовая точка
    convex_shell = [start_dot[2]] #индекс стартовой точки
    ort_base = (x=oneunit(T), y=zero(T)) # начальное направление для поиска следующей точки.

    while next!(convex_shell, dots, ort_base) != start_dot[2]
        ort_base = dots[convex_shell[end]] - dots[convex_shell[end-1]]
    end

	pop!(convex_shell) # удаляем так как будет повторять начальную точку.

    return dots[convex_shell]
end

println("Алгоритм Джарвиса: ", jarvis!( [
		(x=0.0,y=0.0),
		(x=5.0,y=1.0),
		(x=4.0,y=3.0),
		(x=1.0,y=9.0),
		(x=-3.0,y=8.0),
		(x=-5.0,y=2.0),
		(x=-2.0,y=3.0),
	] ) )


# 8. Написать функцию, реализующую алгоритм Грехома построения выпуклой оболочки заданных точек плоскости.
function grekhom!(points::AbstractArray{Vector2D{T}})::AbstractArray{Vector2D{T}} where T
	ydata = [points[i].y for i in 1:length(points)]

	# Сначала надо найти базовую точку , и выбирать базовое направление (точно так же, как это делалось в алгоритме Джарвиса).
	start_dot = findmin(ydata)
	# Далее в выпуклую оболочку помещаются следующие две точки (они гарантированно в неё входят).
	points[begin], points[start_dot[2]] = points[start_dot[2]], points[1]

	# Все остальные точки сортируются по возрастанию угла между вектором и вектором для k = 1,2,3,...,N.
	sort!(@view(points[begin+1:end]), by=(point -> angle(point, (x=oneunit(T),y=zero(T)))))
	push!(points, points[begin])

	convex = [firstindex(points), firstindex(points) + 1, firstindex(points) + 2]

	# Каждая сдедующая в отсортирванном порядке точка помещается в выпуклую оболочку, но пока временно. Т.е. эта точка помещается на вершину стека, в
	# который вконце-концов должна быть помещена вся выпуклая оболочка, но на следующих шагах алгоритма некоторые точки с вершины этого стека могут] быть сняты.
	for i in firstindex(points)+3:lastindex(points)
		while length(convex) > 1 && sign(points[i] - points[convex[end]], points[convex[end-1]] - points[convex[end]]) < 0
			pop!(convex)
		end

		push!(convex, i)
	end

	return points[convex]
end

println("Алгоритм Грехома: ", grekhom!( [
		(x=0.0,y=0.0),
		(x=5.0,y=1.0),
		(x=4.0,y=3.0),
		(x=1.0,y=9.0),
		(x=-3.0,y=8.0),
		(x=-5.0,y=2.0),
		(x=-2.0,y=3.0),
	] ) )


# 9. Написать функцию вычисляющую площадь (ориентированную) заданного многоугольника методом трапеций.

function trapezia_area(poly::AbstractArray{Vector2D{T}})::T where T

    res = zero(T)

	# area = (yk + yk+1)(xk+1 − xk)/2

    for i in firstindex(poly):lastindex(poly)-1
        res .+= (poly[i + 1].y + poly[i]) * (poly[i+1].x - poly[i].x) /2
    end

    return res
end

println("Площадь методом трапеций (квадрата): ", trapezia_area([(x = 0, y = 0), (x = 5, y = 0), (x = 3, y = 5), (x = 0, y = 5)]))

# 10. Написать функцию вычисляющую площадь (ориентированную) заданного многоугольника методом треугольников.

function area_triangle(poly::AbstractArray{Vector2D{T}})::T where T
	
    res = zero(T)

	# area = (yk + yk+1)(xk+1 − xk)/2

    for i in firstindex(poly)+1:lastindex(poly)-1
        res += xdot(poly[i] - (poly[1]), poly[i+1] - poly[1])
    end

    return res
end

println("Площадь методом треугольников: ", area_triangle( [(x = 3.0, y = 1.0), (x = 1.0, y = 2.0), (x = 0.0, y = 1.0), (x = 1.0, y = 0.5)]))

stored_lims = [0,0,0,0]


#Картинка треугольника
draw([(x=0,y=2),
(x=1,y=-1),
(x=-1,y=-1)])

draw((x=0,y=1))

draw((x=2,y=2))

draw((x=-2,y=2))

savefig("C:\\Users\\user\\Desktop\\МОШ\\julia_learning\\tri.png") 
clear()


draw([(x=-5,y=-5),
(x=-5,y=5),
(x=5,y=5),
(x=5,y=-5)])

draw((x=-1,y=0))

draw((x=1,y=0))

savefig("C:\\Users\\user\\Desktop\\МОШ\\julia_learning\\rect.png") 

clear()

poly = [(x=0,y=0)]

for i in 1:15
	p = (x=rand(40:80), y = rand(40:65))

	push!(poly, p)
end

draw(grekhom(poly))
for p in poly
	draw(p)
end
savefig("C:\\Users\\user\\Desktop\\МОШ\\julia_learning\\grekhom.png")
clear()

draw(jarvis!(poly))
for p in poly
	draw(p)
end
savefig("C:\\Users\\user\\Desktop\\МОШ\\julia_learning\\jarvis.png")
clear()