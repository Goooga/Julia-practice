#1
function isprime(n)
    if n<2
        return false
    end
    d = 2
    while d*d <= n
        if mod(n, d) == 0
            return false
        end
        d+=1
    end
    return true
end

#2
function eratosfen(n)
    arr = [i for i in 1:n]
    i = 2
    res=[]
    while i != length(arr)
        if arr[i] != 0
            for j in arr[i]^2:arr[i]:n
                arr[j] = 0
            end
            push!(res, i)
        end
        i+=1
    end
    return res
end

#3
function factor(n)
    arr = eratosfen(Int(n//2))
    res1 = []
    res2 =[]
    for i in arr
        if mod(n, i) == 0
            push!(res1, i)
            k = 0
            while n > 1 && mod(n, i) == 0
                n//=i
                k += 1
            end
            push!(res2, k)
        end
    end
    return res1, res2
end

#4
function meanstd(arr::Union{Tuple, Vector, Array})
    sum1 = 0
    sum2 = 0
    #sum1 - сумма квадратов чисел из массива
    #sum2 = сумма чисел из массива
    for i in arr
        sum1 += i^2
        sum2 += i
    end
    return sqrt(sum1/length(arr) - (sum2/length(arr))^2)
end

# println(meanstd([156, 28, 15, 6, 1, 0]))

#5     
function trace(tree::Vector)
    if isempty(tree)
        return
    end
        
    println(tree[end]) # "обработка" корня
    
    for subtree in tree[1:end-1]
        trace(subtree)
    end
end

#--------------------------------------------------------------
function convert!(intree::Vector, outtree::Dict{Int,Vector})
    #println(outtree)
    if isempty(intree)
        return
    end
    list = []
    for subtree in intree[1:end-1]
        if isempty(subtree)
            push!(list, nothing)
            continue
        end 
        push!(list, subtree[end])
        convert!(subtree,outtree)
    end
    outtree[intree[end]]=list
    return outtree
end

#--------------------------------------

struct Node
    index::Int
    childs::Vector{Node}      
end

function convert(intree::Dict{Int,Vector}, root::Union{Int,Nothing})::Union{Node,Nothing}
    if isnothing(root)
        return nothing
    end
    node = Node(root, [])
    for sub_root in intree[root]
        push!(node.childs, convert(intree, sub_root))
    end
    return node
end
#-------------------------- TEST


tree=Dict{Int, Vector}()

# intree=[[[[],[],6], [], 2], [[[],[],4], [[],[],5], 3],1]
intree = [[[[],[[],[[],[],9],8],7],[[],[],6],3], [[[],[],5],[[],[],4],2], 1]
dict = convert!(intree,tree)
# convert(dict, 1)

function height(tree::Vector)
    if isnothing(tree)
        return 0
    end
    ans=0
    for subtree in tree[1:end-1]
        if height(subtree)>ans
            ans=height(subtree)
        end
    end
    return 1 + ans
end

# println(convert!(intree, dict))
# println(height(intree))

function count_leafs(tree::Vector)
    if length(tree) == 0
        return 0
    end

    ans=0
    isleaf = true
    for subtree in tree[1:end-1]
        if length(subtree) != 0
            isleaf = false
            break
        end
    end

    if isleaf == true
        return 1
    else
        for subtree in tree[1:end - 1]
            ans += count_leafs(subtree)
        end
    end
    return ans
end

# println(count_leafs(intree))
function nodes_count(intree)
    node_count = 1

    is_leaf = true
    for subtree in intree[1:end-1]
        is_leaf = false
        node_count += nodes_count(subtree)
    end
    if length(intree) == 0
        return 0
    end
    if is_leaf
        return 1
    end
    
    return node_count

end

function max_valence(intree)
    valence = 0
    arr_valence = []
    for subtree in intree[1:end-1]
        valence += 1
        push!(arr_valence, (max_valence(subtree)))
    end
    push!(arr_valence, valence)
    return maximum(arr_valence)
end

function sum_ways(intree)
    sum_ways = 0
    for subtree in intree[1:end-1]
        sum_ways+= 1
        sum_ways += sum_ways(subtree)
    end
    return sum_ways
end

function average_lenght(intree)
    N = nodes_count(intree)
    S = sum_ways(intree)
    return (S + N)/N
end

# println(max_valence(intree))