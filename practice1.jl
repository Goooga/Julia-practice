import Base.display

#1 НОД
function gcd_(a::T, b::T) where T
    while !iszero(b)
        a, b = b, mod(a, b)
    end
    return a
end

function remdiv(a::T, b::T) where T
    return rem(a, b), div(a,b)
end

#2 gcdx_
function gcdx_(a::T, b::T) where T
   u, v = one(T), zero(T)
   u′, v′ = v, u
    # ИНВ-т d = ua + vb 
   while !iszero(b)
        n, m = remdiv(a,b)
        a, b = b, n
        u, v = v, u - m*v
        u′, v′ = v′, v - m*v′
   end
   if isnegative(a)
     a, u, v = -a, -u, -v
   end
   return a, u, v
end

# function gcdx_(a::Int, b::Int)
#     a, u, v = gcdx_(a, b)
#     if a < 0
#         a, u, v = -a, -u, -v
#     end
# end

isnegative(a::Number) = (a<0)

#3 обратный в кольце
function invmod_(a::T, M::T) where T
    if gcd(a, M) != 1
        return nothing
    end
    res = gcdx_(a, M)
    return mod(res[2],M)
end

#4 диафантово уравнение
function diaphant_solve(a::T,b::T,c::T) where T
    #xa + yb = 1 
    
    gcd = gcd_(a, b)
    while gcd != 1
        if gcd > 1 && mod(c, gcd) != 0
            return nothing
        end
        a, b, c = div(a,gcd), div(b,gcd), div(c,gcd)
        gcd = gcd_(a, b)
    end
    return gcdx_(a, b)[2] * c, gcdx_(a, b)[3] * c
end

#5 кольцо вычетов
struct Residue{T,M}
    a::T
    Residue{T,M}(a::T) where{T,M} = new(mod(a, M))
    #function Residue{Polynom, M}(a::Polynom) where{M} 
end

Base.zero(::Type{Residue{T, M}}) where {T, M} = Residue{T, M}(zero(T)) 

Base.zero(a::Residue{T,M}) where{T,M} = Residue{T,M}(zero(T))

function Base.:+(x::Residue{T, M}, y::Residue{T, M}) where{T,M}
    return Residue{T, M}(x.a + y.a)
end

function  Base.:+(x::Residue{T, M}, y::Int) where {T, M}
    return Residue{T, M}(x.a + y)
end

function Base.:-(x::Residue{T, M}) where {T, M}
    return Residue{T, M}(-(x.a))
end

function Base.:-(x::Residue{T, M}, y::Residue{T, M}) where {T,M}
    return x + (-y)
end

function Base.:-(x::Residue{T,M}, y::Int) where{T,M}
    return Residue{T,M}(x.a - y)
end

function Base.:^(x::Residue{T,M}, y::Int) where{T,M}
    return Residue{T,M}(x.a^y)
end

function Base.:^(x::Residue{T,M}, y::Residue{T,M}) where{T,M}
    return Residue{T,M}(x.a^y.a)
end

function Base.:*(x::Residue{T,M}, y::Int) where{T,M}
    return Residue{T,M}(x.a*y)
end

function Base.:*(x::Residue{T,M}, y::Residue{T,M}) where{T,M}
    return Residue{T,M}(x.a*y.a)
end

function inverse(x::Residue{T,M}) where {T,M}
    if (diaphant_solve(x.a,-M,1)) != nothing
        k1,_=diaphant_solve(x.a,-M,1)
        return Residue{T,M}(k1)
    end
    return nothing
end

Base.zero(::Type{Residue{T, M}}) where {T, M} = zero(T)

function Base.display(x::Residue{T,M}) where {T,M}
    return(display(x.a))
end


#6 Полином
struct Polynom{T,Deg} 
    cfs::AbstractVector{T}
    Polynom{T,Deg}() where {T,Deg} = new(Vector{T}(undef,Deg+1) )
end

function make_polynom(a::Union{Vector{T}, Tuple{T}}) where{T}
    len=length(a)
    while (a[len]==zero(T)) && (len>1)
        len-=1
    end
    res=Polynom{T,len-1}()
    for i in 1:len
        res.cfs[i]=a[i]
    end
    return res
end

get_degree(x::Polynom{T, Deg}) where {T, Deg} = Deg

function Base.:+(a::Polynom{T, Deg1}, b::Polynom{T, Deg2}) where {T, Deg1, Deg2}
    if Deg1 < Deg2
        a, b = b, a
    end
    arr = Vector{T}(undef, max(Deg1, Deg2)+1)
    for i in 1:max(Deg1, Deg2) + 1
        arr[i] = a.cfs[i]
    end
    for i in 1:min(Deg1, Deg2) + 1
        arr[i]+=b.cfs[i]
    end
    return make_polynom(arr)
end

function Base.:-(a::Polynom{T, Deg}) where {T, Deg}
    res = Polynom{T, Deg}()
    for i in 1:Deg+1
        res.cfs[i] = -a.cfs[i]
    end
    return res
end

function Base.:-(a::Polynom{T, Deg1}, b::Polynom{T, Deg2}) where {T, Deg1, Deg2}
    return a +(- b)
end

function Base. *(x::Polynom{T,D1}, y::Polynom{T,D2}) where {T,D1,D2}
    res=Vector{T}(undef,D1+D2+1)
    for i in 0:D1+D2
        c=0
        for j in 0:i
            if 0<j+1<=D1+1 && 0<i-j+1<=D2+1
                c+=x.cfs[j+1]*y.cfs[i-j+1]
            end
        end
        res[i+1]=c
    end
    return make_polynom(res)
end

function Base. mod(a::Polynom{T,D1}, b::Polynom{T,D2}) where {T,D1,D2}

    while get_degree(a)>=get_degree(b)
        x=Polynom{T,get_degree(a)-get_degree(b)}()
        x.cfs[get_degree(a)-get_degree(b)+1]=a.cfs[get_degree(a)+1]/b.cfs[get_degree(b)+1]
        a-=x*b

    end
    return a
end

function Base. mod(x::Polynom{T,D1}, t::Tuple) where {T,D1}
    v=[i for i in t]
    y=make_polynom(v)
    while get_degree(x)>=get_degree(y)
        a=Polynom{T,get_degree(x)-get_degree(y)}()
        a.cfs[get_degree(x)-get_degree(y)+1]=x.cfs[get_degree(x)+1]/y.cfs[get_degree(y)+1]

        x-=y*a

    end
    return x
end


function Base.display(a::Polynom{T, Deg}) where{T, Deg}
    res=""
    for i in 0:Deg
        res*= string(a.cfs[i + 1])*" t^"*string(i)*" + "
    end
    return res
end

function Base.one(::Type{Polynom{T}}) where {T}
    res = make_polynom(Tuple(one(T)))
    return res
end

# function Base.abs(p::Polynom{T, Deg}) where {T, Deg}
#     res = make_polynom(Vector{T}(undef, Deg + 1))
#     for i in 1:Deg + 1
#         res.cfs[i] = abs(p.cfs[i])
#     end
#     return res
# end

function polyval(p::Polynom, x)
    n = length(p.cfs)
    res = p.cfs[n]
    res_d = 0
    while n != 1
        #println(n)
        n -= 1
        res_d = res_d * x + res
        res *= x
        res += p.cfs[n]# res'(x) = res'(x)*x + res(x)
        #println(res)
    end
    return res, res_d
end

(p::Polynom)(x) = polyval(p, x)


#display(one(Polynom{Int}))

# print(gcdx_(-13, -5))

# p = make_polynom([9, 7, 5, 3, 1])
# p(2)

#8 взаимодействие полинома и кольца
# C=make_polynom([8,2.,4.,5.])
 #D=make_polynom([9.,7.,15.])
# a=Residue{Polynom,(0.,0.,1.)}(C)
# println(mod(C,(0.,0.,1.)))
# println(display(a))

# arr = [Residue{Int, 3}(6), Residue{Int, 3}(-8), Residue{Int, 3}(2)]
# p = make_polynom(arr)
# print(display(p))
