
# № 1 Размещения n элементов по k с повторениями.

function next_repit_placement!(p::Vector{T}, n::T) where T<:Integer
    i = findlast(x->(x < n), p) # используется встроенная функция высшего порядка
    # i - это последний(первый с конца) индекс: x[i] < n, или - nothing, если такого индекса нет (p == [n,n,...,n])
    isnothing(i) && (return nothing)
    p[i] += 1
    p[i+1:end] .= 1 # - устанавливаются минимально-возможные значения
    return p
end

#=Данный код определяет функцию next_repit_placement!, 
которая принимает вектор p и целочисленный аргумент n. Функция изменяет вектор p следующим образом:
находит последний (первый с конца) индекс i, для которого p[i] < n,
увеличивает значение p[i] на 1 и устанавливает все последующие элементы вектора p равными 1.
Если такого индекса i не существует (вектор p состоит только из значений n), то функция возвращает nothing.=#


# println("Генерация всех размещений с повторениями из n элементов {1,2,...,n} по k")
# p = [1, 1, 1]
# print(p)
function test(n=2, k=3)
    p = ones(Int, k)
    while !isnothing(p)
        println(p)
        p = next_repit_placement!(p, 3)
    end
end


# № 2 Перестановки n элементов
function next_permute!(p::AbstractVector)
    n = length(p)
    k = 0 # firstindex(p) - 1
    for i in reverse(1:n - 1) # reverse(firstindex(p):lastindex(p) - 1)
        if p[i] < p[i + 1]
            k = i
            break
        end
    end
    k == firstindex(p) - 1 && return nothing # p[begin] > p[begin + 1] > ... > p[end]

    #утв: p[k] < p[k + 1] > p[k + 2] > ... > p[end]
    i = k + 1
    while i < n && p[i + 1] > p[k] # i < lastindex(p) && p[i + 1] > p[k]
        i += 1
    end
    #утв: p[i] - наименьшее из p[k + 1:end], большее p[k]
    p[k], p[i] = p[i], p[k]
    #утв: по-прежнему p[k + 1]>...>p[end]
    reverse!(@view p[k + 1:end])
    return p
end
#=
Алгоритм работы функции заключается в следующем:
1. Находим индекс k такой, что p[k] < p[k+1], если такого индекса нет, то перестановка последняя в лексикографическом порядке.
2. Находим индекс i такой, что p[i] - наименьшее из p[k+1:end], большее p[k].
3. Меняем местами p[k] и p[i].
4. Разворачиваем отрезок p[k+1:end].
=#
# println("Генерация всех размещений с повторениями из n элементов {1,2,...,n} по k")
p = [1, 2, 3]
# while !isnothing(p)
#     println(p)
#     p = next_permute!(p)
# end

# № 3 Все подмножества n -элементного множества
# println("Генерация всех всех подмножеств n-элементного множества {1,2,...,n} 1 способ")
# Первый способ - на основе генерации двоичных кодов чисел 0, 1, ..., 2^n-1

indicator(i::Integer, n::Integer) = reverse(digits(Bool, i; base=2, pad=n))

# println("1 способ")
# println(indicator(12, 5))

# Второй способ - на основе непосредственной генерации последовательности индикаторов в лексикографическом порядке

function next_indicator!(indicator::AbstractVector{Bool})
    i = findlast(x->(x==0), indicator)
    isnothing(i) && return nothing
    indicator[i] = 1
    indicator[i+1:end] .= 0
    return indicator 
end

# println("2 способ")
# println(next_indicator!(indicator(12, 5)))


# № 4 k -элементные подмножества n -элементного множества.

function next_indicator!(indicator::AbstractVector{Bool}, k)
    # в indicator - ровно k единц, остальные - нули, но это не проверяется! (фактически k - не используется)
    i = lastindex(indicator)
    while indicator[i] == 0
        i -= 1
    end
    #УТВ: indic[i] == 1 и все справа - нули
    m = 0 
    while i >= firstindex(indicator) && indicator[i] == 1 
        m += 1
        i -= 1
    end
    if i < firstindex(indicator)
        return nothing
    end
    #УТВ: indicator[i] == 0 и справа m > 0 единиц, причем indicator[i + 1] == 1
    indicator[i] = 1
    indicator[i + 1:end] .= 0
    indicator[lastindex(indicator) - m + 2:end] .= 1
    return indicator 
end

# println("Генерация всех k-элементных подмножеств n-элементного множества {1, 2, ..., n}")
n = 6
k = 3
a= 1:6

# a[findall(next_indicator!([zeros(Bool, n-k); ones(Bool, k)], k))] |> println


# № 5 Разбиения натурального числа на положительные слагаемые.

function next_split!(s ::AbstractVector{Int64}, k)
    k == 1 && return (nothing, 0)
    i = k-1 # - это потому что s[k] увеличивать нельзя
    while i > 1 && s[i-1]>=s[i]
        i -= 1
    end
    #УТВ: i == 1 или i - это наименьший индекс: s[i-1] > s[i] и i < k
    s[i] += 1
    #Теперь требуется s[i+1]... - уменьшить минимально-возможным способом (в лексикографическом смысле) 
    r = sum(@view(s[i+1:k]))
    k = i+r-1 # - это с учетом s[i] += 1
    s[(i+1):(length(s)-k)] .= 1
    return s, k
end

# println("Генерация всех разбиений натурального числа на положительные слагаемые")
# println(next_split!(ones(Int64, 5), 5))


# Абстрактный пользовательский тип для генерации комбинаторных объектов
abstract type AbstractCombinObject
    # value::Vector{Int} - это поле предполагается у всех конкретных типов, наследующих от данного типа
end

Base.iterate(obj::AbstractCombinObject) = (get(obj), nothing)
Base.iterate(obj::AbstractCombinObject, state) = (isnothing(next!(obj)) ? nothing : (get(obj), nothing))

# № 6 Итераторы для перебора комбинаторных объектов.
# Размещения с повторениями
struct RepitPlacement{N,K} <: AbstractCombinObject
    value::Vector{Int}
    RepitPlacement{N,K}() where {N, K} = new(ones(Int, K))
end

Base.get(p::RepitPlacement) = p.value
next!(p::RepitPlacement{N,K}) where {N, K} = next_repit_placement!(p.value, N)

# Тестирование:
println("Размещения с повторениями")
for a in RepitPlacement{2,3}() 
    println(a)
end


#Перестановки
struct Permute{N} <: AbstractCombinObject
    value::Vector{Int}
    Permute{N}() where N = new(collect(1:N))
end

Base.get(obj::Permute) = obj.value
next!(permute::Permute) = next_permute!(permute.value)

# Тест:
println("Перестановки")
for p in Permute{4}()
    println(p)
end

#Все подмножества N-элементного множества
struct Subsets{N} <: AbstractCombinObject
    indicator::Vector{Bool}
    Subsets{N}() where N = new(zeros(Bool, N))
end

Base.get(sub::Subsets) = sub.indicator
next!(sub::Subsets) = next_indicator!(sub.indicator) 

#Тест:
println("Все подмножества N-элементного множества")
for sub in Subsets{3}()
    println(sub)
end

#k-элементные подмоножества n-элементного множества
struct KSubsets{M,K} <: AbstractCombinObject
    indicator::Vector{Bool}
    KSubsets{M, K}() where{M, K} = new([zeros(Bool, length(M)-K); ones(Bool, K)])
end

Base.get(sub::KSubsets) = sub.indicator
next!(sub::KSubsets{M, K}) where{M, K} = next_indicator!(sub.indicator, K) 

println("3-элементные подмножества множества из 4 элементов")
#Тест:
for sub in KSubsets{1:4, 3}()
    sub |> println
end

#Разбиения
mutable struct NSplit{N} <: AbstractCombinObject
    value::Vector{Int}
    num_terms::Int # число слагаемых (это число мы обозначали - k)
    NSplit{N}() where N = new(vec(ones(Int, N)), N)
end

Base.get(nsplit::NSplit) = nsplit.value[begin:nsplit.num_terms]
function next!(nsplit::NSplit)
    a, b = next_split!(nsplit.value, nsplit.num_terms)
    if isnothing(a) return nothing end
    nsplit.value, nsplit.num_terms = a, b
    get(nsplit)
end

# Тест:
println("Разбиения 5")
for s in NSplit{5}()
    println(s)
end


function dfs(graph::AbstractDict, vstart::T) where T <: Integer
    stack = [vstart]
    mark = zeros(Bool, length(graph)) # length(graph) = N = число всех вершин графа
    mark[vstart] = true
    while !isempty(stack)
        v = pop!(stack)
        for u in graph[v] mark[u] = (!mark[u] ? (push!(stack, u); true) : true) end
    end
    return mark
end


graph1 = Dict{Int64, Vector{Int64}}([(1, [3]), (2, [4]), (3, [1]), (4, [2, 5]), (5, [4])])
graph2 = Dict{Int64, Vector{Int64}}([(1, [2, 3]), (2, [1, 4]), (3, [1]), (4, [2, 5]), (5, [4])])
graph3 = Dict{Int64, Vector{Int64}}([(1, [2, 3]), (2, [1, 4]), (3, [1, 6]), (4, [2, 5]), (5, [4, 6]), (6, [3, 5])])
graph4 = Dict{Int64, Vector{Int64}}([(1, [2, 3]), (2, [1, 3, 4]), (3, [1, 2]), (4, [2, 5]), (5, [4])])
println(dfs(graph1, 1))

function isConnectedGraph(graph::AbstractDict) :: Bool
    res = dfs(graph, 1)
    return (sum(res) == length(res) ? true : false)
end

println(isConnectedGraph(graph1))

function find_connectivity_components(graph::AbstractDict, len = length(graph))
    mark = ones(Bool, len)
    ans = []
    for i in 1:len
        res = mark[i] ? (t = dfs(graph, i);  push!(ans, t); t) : Bool[0]
        mark[findall(res)] .= false
    end
    return ans
end

println(find_connectivity_components(graph1))


function  isBipartite(graph::AbstractDict, len = length(graph)) :: Bool
    color = fill(-1, len)
    queue = []
    for i in 1:len
        if color[i] != -1 continue end
        color[i] = 0
        push!(queue, i)
        while !isempty(queue)
            v = popfirst!(queue)
            if !isnothing(findfirst(isequal(color[v]), color[graph[v]])) return false end
            found = graph[v][findall(isequal(-1), color[graph[v]])]
            color[found] .= (color[v] + 1) % 2
            append!(queue, found)
        end
    end
    return true
end


println(isBipartite(graph3))
println(isBipartite(graph4))

