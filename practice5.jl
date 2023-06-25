import Base.sort!
import Base.sort

function bubble_sort!(a::AbstractArray)
    n = length(a)
    for i in 1:n
        for j in i:n
            if a[i]>a[j]
                a[i], a[j] = a[j], a[i]
                
            end
        end
    end
    return a
end

bubble_sort(b::AbstractArray)=bubble_sort!(copy(b))
    

function comb_sort!(a::AbstractArray; factor=1.2473309)
    step = length(a)
    while step >= 1
        for i in 1:length(a)-step
            if a[i] > a[i+step]
                a[i], a[i+step] = a[i+step], a[i]
            end
        end
        step = Int(floor(step/factor))
    end
    bubble_sort!(a)
end
comb_sort(b::AbstractArray; factor=1.2473309)=comb_sort!(copy(b); factor=factor)
    
function part_sort!(A, b)
    N = length(A)
    K=0
    L=0
    M=N
   
    while L < M 
        if A[L+1] == b
            L += 1
        elseif A[L+1] > b
            A[L+1], A[M] = A[M], A[L+1]
            M -= 1
        else 
            L += 1; K += 1
            A[L], A[K] = A[K], A[L]
        end
    end
    return K, M+1 

end
part_sort(A, b)=part_sort!(copy(A, b))
# ∀i∈1:K   A[i]<b
# \forall i \in K+1:L \ \ \ A[i]==b∀i∈K+1:L   A[i]==b
# \forall i \in M+1:N \ \ \ A[i]>b∀i∈M+1:N   A[i]>b

function quick_sort!(A)
    if isempty(A)
        return A
    end
    N = length(A)
    K, M = part_sort!(A, A[rand(1:N)]) 
    quick_sort!(@view A[1:K])
    quick_sort!(@view A[M:N])
    return A
end

quick_sort(A)=quick_sort!(copy(A))

function sort_perm!(a)
    indexes = collect(firstindex(a):lastindex(a))
    n = length(a)
    for i in 1:n
        for j in i:n
            if a[i]>a[j]
                a[i], a[j] = a[j], a[i]
                indexes[i], indexes[j] = indexes[j], indexes[i]
            end
        end
    end
    return indexes
end

sort_perm(a)=sort_perm!(copy(a))




a = [1, 4, 6, 2, 3, 5]

# println(sort_perm(a))
# println(a)

function shell_sort!(a; step_series = (length(a)÷2^i for i in 1:Int(floor(log2(length(a))))) )
    for step in step_series
        for i in firstindex(a):lastindex(a)-step
            j = i
            while j >= firstindex(a) && a[j] > a[j+step]
                a[j], a[j+step] = a[j+step], a[j]
                j -= step
            end
        end
    end
    return a
end

#4
function insert_sort!(vector)
    n = 1
    while n < length(vector) 
        n += 1
        i = n
        while i > 1 && vector[i-1] > vector[i]
            vector[i], vector[i-1] = vector[i-1], vector[i]
            i -= 1
        end
    
    end
    return vector
end

insert_sort(vector)=insert_sort!(copy(vector))



sort!(lst)=insert_sort!(lst)

sort(lst)=sort!(copy(lst))


#5
function calc_sort!(A::AbstractVector{<:Integer})
    min_val, max_val = extrema(A)
    num_val = zeros(Int, max_val-min_val+1) # - число всех возможных значений
    for val in A
        num_val[val-min_val+1] += 1
    end  
    k = 0
    for (i, num) in enumerate(num_val)
        A[k+1:k+num] .= min_val+i-1
        k += num
    end
    return A
end

n = 100000
lst = randn(n)
a = copy(lst)
b = copy(lst)
println("comb_sort! time:")
@time begin
    comb_sort!(a)

end
@time begin
    shell_sort!(b)
end
c = copy(lst)
println("insert_sort time:")
@time begin
    insert_sort!(c)
end


@inline function Base.merge!(a1, a2, a3)::Nothing # @inline - делает функцию "встраиваемой", т.е. во время компиляции ее тело будет встроено непосредственно в код вызывающей функции (за счет этого происходит экономия на времени, затрачиваемым на вызов функции; это время очень небольшое, но тем не менее)
    i1, i2, i3 = 1, 1, 1
    @inbounds while i1 <= length(a1) && i2 <= length(a2) # @inbounds - передотвращает проверки выхода за пределы массивов
        if a1[i1] < a2[i2]
            a3[i3] = a1[i1]
            i1 += 1
        else
            a3[i3] = a2[i2]
            i2 += 1
        end
        i3 += 1
    end
    @inbounds if i1 > length(a1)
        a3[i3:end] .= @view(a2[i2:end]) # Если бы тут было: a3[i3:end] = @view(a2[i2:end]), то это привело бы к лишним аллокациям (к созданию промежуточного массива)
    else
        a3[i3:end] .= @view(a1[i1:end])
    end
    nothing
end
function merge_sort!(a)
    b = similar(a) # - вспомогательный массив того же размера и типа, что и массив a
    N = length(a)
    n = 1 # n - текущая длина блоков
    @inbounds while n < N
        K = div(N,2n) # - число имеющихся пар блоков длины n
        for k in 0:K-1
            merge!(@view(a[(1:n).+k*2n]), @view(a[(n+1:2n).+k*2n]), @view(b[(1:2n).+k*2n]))
        end
        if N - K*2n > n # - осталось еще смержить блок длины n и более короткий остаток
            merge!(@view(a[(1:n).+K*2n]), @view(a[K*2n+n+1:end]), @view(b[K*2n+1:end]))
        elseif 0 < N - K*2n <= n # - оставшуюся короткую часть мержить не с чем
            b[K*2n+1:end] .= @view(a[K*2n+1:end])
        end
        a, b = b, a
        n *= 2
    end
    if isodd(log2(n)) # - если цикл был выполнен нечетное число раз, то b - это исходная ссылка на массив (на внешний массив), и a - это ссылка на вспомогательный массив (локальный)
        b .= a # b = copy(a) - это было бы не то же самое, т.к. при этом получилась бы ссылка на новый массив, который создает функция copy
        a = b
    end
    return a # - исходная ссылка на внешний массив (проверить, что это так, можно с помощью ===)
end
d = copy(lst)
println("merge time:")
@time begin
    merge_sort!(d)
end

function part_sort!(A, b)
    N = length(A)
    K=0
    L=0
    M=N
    #ИНВАРИАНТ: A[1:K] < b && A[K+1:L] == b && A[M+1:N] > b
    while L < M 
        if A[L+1] == b
            L += 1
        elseif A[L+1] > b
            A[L+1], A[M] = A[M], A[L+1]
            M -= 1
        else # if A[L+1] < b
            L += 1; K += 1
            A[L], A[K] = A[K], A[L]
        end
    end
    return K, M+1 
    # 1:K и M+1:N - эти диапазоны индексов определяют ещё не 
    # отсортированные части массива A
end
function findMedian(a, n)
    if (n % 2 == 0)
        part_sort!(a, n÷2)
        part_sort!(a, (n - 1) ÷ 2)
        return (a[(n - 1) ÷ 2] + a[n ÷ 2]) ÷ 2
    else 
        part_sort!(a, n÷2)
        return a[n ÷ 2]
    end
end
f = rand(Int, 1000)
println("median time:")
@time findMedian(f, length(f))

N = Int64(1e7)
a = rand(1:10000,N)

function linear_sort(a::AbstractVector{T}, N::Int) where T
    arr = zeros(N)
    res = []
    for i in 1:n
        arr[a[i]] += 1
    end
    for i in 1:n
        while arr[i]>0
            push!(res,i)
            arr[i]-=1
        end
    end
    return c
end

println("linear_sort time:")
@time linear_sort(a, N)