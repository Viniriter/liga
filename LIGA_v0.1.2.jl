import Base.show, Base.copy, Base.subtypes
importall Base.Operators
#########################################################
#########################################################
"""
```
kb(x::Vector{Bool}, a::Number)
```
Creates a scaled basis element of geometric algebra space from x and a.

```julia-repl
julia> kb(e13, 2.3)
2.3e13
julia> kb(e23)
1.0e23
```
"""
type kb
    e::Vector{Bool}
	scl::Number

	kb(a::Vector{Bool}) =   new(a, 1.0)
	kb(a::Vector{Bool}, b::Number) =   new(a, b)
end
#########################################################
"""
```
kMultvec(X::Vector{kb})
```

Creates a multivector in the Euclidean geometric algebra (Gₚ) composed by the sum of the coordinates of X.

```julia-repl
julia> kMultvec([kb(e1, 2.0), kb(e12, 3.3), kb(e123)])
2.0e1 + 3.3e12 + 1.0e123
```
"""
type kMultvec
    comp::Vector{kb}

	kMultvec(a::kb) = new([a])
	kMultvec(A::Vector{kb}) = new(A)
	kMultvec(a::Number) = new([kb(id, a)])
end
#########################################################
"""
```
grade(a::kb)
grade(a::pb)
grade(a::cb)
```
Returns the grade of the basis element x

```julia-repl
julia> grade(kb(e13, 1.5))
2
```
"""
function grade(a::kb)
	acu = 0
	for i=1:length(a.e)
		if a.e[i] == true
			acu = acu + 1
		end
	end
	return acu
end
#########################################################
"""
```
mvectovec(X::Vector{kMultvec})
```
Auxiliary function.

"""
function mvectovec(X::Vector{kMultvec})
	l = length(X)
	flag = 0
	for i=1:l
		m = length(X[i].comp)
		for j=1:m
			if grade(X[i].comp[j]) != 1
				flag = 1
			end
		end
	end
	if flag == 1
		error("Grade error!")
		return 0
	else
		n = length(X[1].comp[1].e)
		A = Matrix{Number}(l, n)
		fill!(A, 0.0)
		for i=1:l
			m = length(X[i].comp)
			for j=1:m
				for k=1:n
					if X[i].comp[j].e[k] == true
						A[i,k] = X[i].comp[j].scl
					end
				end
			end
		end
		return A
	end
end
#########################################################
"""
```
kblade(X::Vector{kMultvec})
```
Creates a blade in the Euclidean geometric algebra (Gp) composed by the outer product of the elements of X.

The blade is created only if the multivectors that compound X are both 1-vectors and X is a L.i set

```julia-repl
julia> a = kMultvec([kb(e1,2.0),kb(e3,4.0),kb(e3,4.0)])
2.0e1 + 4.0e3 + 4.0e3
julia> b = kMultvec([kb(e1, 1.0), kb(e2, 1.0), kb(e3,2.0)])
1.0e1 + 1.0e2 + 2.0e3
julia> A = kblade([a,b])
(2.0e1 + 4.0e3 + 4.0e3)∧(1.0e1 + 1.0e2 + 2.0e3)
```

"""
type kblade
	conj::Vector{kMultvec}

	function kblade(X::Vector{kMultvec})
		T = mvectovec(X)
		if det(T*transpose(T)) != 0
			new(X)
		else
			error("L.d")
		end
	end
end
#########################################################
function copy(a::kb)
	e = copy(a.e)
	esc = copy(a.scl)
	result = kb(e, esc)
	return result
end
#########################################################
"""
```
geoprod(a::kb, b::kb)
geoprod(X::kMultvec, Y::kMultvec)
geoprod(A::kblade, B::kblade)
geoprod(a::pb, b::pb)
geoprod(X::pMultvec, Y::pMultvec)
geoprod(A::pblade, B::pblade)
geoprod(a::cb, b::cb)
geoprod(X::cMultvec, Y::cMultvec)
geoprod(A::cblade, B::cblade)

```
Receives as input parameters elements of the same type and returns the geometric product between them.

It can be used the operator "∘ (circ)" instead of call the function.
"""
function geoprod(a::kb, b::kb)
	l = length(a.e)
	Ea = copy(a.e)
	Eb = copy(b.e)
	ab = kb(copy(Ea), 0.0)
	scl = 1
	eq = 0
	cont = 0
	for i=1:l-1
		if Eb[i] == true
			for j=i+1:l
				if Ea[j] == true
					cont = cont + 1
				end
			end
		end
	end
	for i=1:l
		if Ea[i] != Eb[i]
			ab.e[i] = true
		else
			ab.e[i] = false
		end
	end
	ab.scl = a.scl*b.scl*(-1)^cont
	return ab
end
#########################################################
function Base.:∘(a::kb,b::kb)
    return geoprod(a,b)
end
#########################################################
"""
```
inner(a::kb, b::kb)
inner(X::kMultvec, Y::kMultvec)
inner(A::kblade, B::kblade)
inner(a::pb, b::pb)
inner(X::pMultvec, Y::pMultvec)
inner(A::pblade, B::pblade)
inner(a::cb, b::cb)
inner(X::cMultvec, Y::cMultvec)
inner(A::cblade, B::cblade)

```
Receives as input parameters elements of the same type and returns the inner product between them.

It can be used the operator "⋅ (cdot)" instead of call the function.
"""
function inner(a::kb, b::kb)
	if grade(geoprod(a,b)) == abs(grade(a) - grade(b))
		return geoprod(a,b)
	else
		x = Vector{Bool}(length(a.e))
		fill!(x, false)
		return kb(x, 0.0)
	end
end
#########################################################
function Base.:⋅(a::kb,b::kb)
    return inner(a,b)
end
#########################################################
"""
```
outer(a::kb, b::kb)
outer(X::kMultvec, Y::kMultvec)
outer(A::kblade, B::kblade)
outer(a::pb, b::pb)
outer(X::pMultvec, Y::pMultvec)
outer(A::pblade, B::pblade)
outer(a::cb, b::cb)
outer(X::cMultvec, Y::cMultvec)
outer(A::cblade, B::cblade)

```
Receives as input parameters elements of the same type and returns the outer product between them.

It can be used the operator " ^ " instead of call the function.
"""
function outer(a::kb, b::kb)
	if grade(geoprod(a,b)) == (grade(a) + grade(b))
		return geoprod(a,b)
	else
		x = Vector{Bool}(length(a.e))
		fill!(x, false)
		return kb(x, 0.0)
	end
end
#########################################################
function Base.:^(a::kb,b::kb)
    return outer(a,b)
end
#########################################################
"""
```
scalar(a::kb, b::kb)
scalar(X::kMultvec, Y::kMultvec)
scalar(A::kblade, B::kblade)
scalar(a::pb, b::pb)
scalar(X::pMultvec, Y::pMultvec)
scalar(A::pblade, B::pblade)
scalar(a::cb, b::cb)
scalar(X::cMultvec, Y::cMultvec)
scalar(A::cblade, B::cblade)

```
Receives as input parameters elements of the same type and returns the scalar product between them.

It can be used the operator " * " instead of call the function.
"""
function scalar(a::kb, b::kb)
	if grade(geoprod(a,b)) == 0
		return geoprod(a,b).scl
	else
		x = Vector{Bool}(length(a.e))
		fill!(x, false)
		return kb(x, 0.0)
	end
end
#########################################################
function Base.:*(a::kb,b::kb)
    return scalar(a,b)
end
#########################################################
function bscalar(a::kb, b::kb)
	if grade(geoprod(a,b)) == 0
		return geoprod(a,b)
	else
		x = Vector{Bool}(length(a.e))
		fill!(x, false)
		return kb(x, 0.0)
	end
end
#########################################################
"""
```
mvsum(a::kb, b::kb)
mvsum(X::kMulvec, Y::kMultvec)
mvsum(a::pb, b::pb)
mvsum(X::pMulvec, Y::pMultvec)
mvsum(a::cb, b::cb)
mvsum(X::cMulvec, Y::cMultvec)


```
Receives as input parameters elements of the same type and returns the sum of them.

It can be used the operators "+" or "-" instead of call the function.
"""

function mvsum(a::kb, b::kb)
	l = length(a.e)
	e = Vector{Bool}(l)
	fill!(e, false)
	sume = kb(e, 0.0)
	if a.e == b.e
		sume.e = a.e
		sume.scl = a.scl + b.scl
	elseif a.e != b.e
		sume = kMultvec([a,b])
	end
	return sume
end
#########################################################
function Base.:+(a::kb,b::kb)
    return mvsum(a,b)
end
function Base.:-(a::kb,b::kb)
	b.scl = -b.scl
	return mvsum(a,b)
end
#########################################################
"""
```
mvreverse(a::kb)
mvreverse(X::kMultvec)
mvreverse(A::kblade)
mvreverse(a::pb)
mvreverse(A::pblade)
mvreverse(a::cb)
mvreverse(A::cblade)

```
Returns the reverse of the element.
"""
function mvreverse(a::kb)
	k = grade(a)
	e = copy(a.e)
	b = kb(e)
	b.scl = (-1)^(k*(k-1)/2)*a.scl
	return b
end
#########################################################
"""
```
dual(a::kb)
dual(X::kMultvec)
dual(A::kblade)
dual(X::pMultvec)
dual(a::pb)
dual(A::pblade)

```
Returns the dual element of the basis blade a, the multivetor X or the blade A.
"""
function dual(a::kb)
	l = length(a.e)
	e = Vector{Bool}(l)
	fill!(e, true)
	I = kb(e)
	Inv = mvreverse(I)
	dual = geoprod(a, Inv)
	return dual
end
#########################################################
"""
```
magnitude(a::kb)
magnitude(a::kMultvec)
magnitude(A::kblade)
magnitude(a::pb)
magnitude(a::pMultvec)
magnitude(A::pblade)
magnitude(a::cb)
magnitude(A::cblade)


```
Returns the magnitude of the input element, that is.
"""
function magnitude(a::kb)
	b = mvreverse(a)
	return sqrt(scalar(a,b))
end
#########################################################
"""
```
remove(A::kMultvec, b::Number)
remove(A::pMultvec, b::Number)
remove(A::cMultvec, b::Number)

```
Auxiliary function.
"""
function remove(A::kMultvec, b::Number)
	l = length(A.comp)
	aux = 0
	for i=1:l
		if A.comp[i].scl == b
			aux = aux+1
		end
	end
	B = Vector{kb}(l-aux)
	for i=1:l
		for j=1:l-1
			if A.comp[j].scl == b
				aux2 = A.comp[j+1]
				A.comp[j+1] = A.comp[j]
				A.comp[j] = aux2
			end
		end
	end
	for i=1:(l-aux)
		B[i] = A.comp[i]
	end
	return kMultvec(B)
end
#########################################################
function mvsum(A::kMultvec, B::kMultvec)
	l = length(A.comp)
	m = length(B.comp)
	AB = Vector{kb}(l+m)
	for i=1:l
		AB[i] = copy(A.comp[i])
	end
	for i=l+1:l+m
		AB[i] = copy(B.comp[i-l])
	end
	for j=1:(m+l-1)
		aux = AB[j]
		for i=j+1:m+l
			if AB[i].e == aux.e && AB[i].scl != 0.0
				AB[j] = mvsum(aux, AB[i])
				aux = AB[j]
				AB[i].scl = 0.0
			end
		end
	end
	result = remove(kMultvec(AB), 0.0)
	if length(result.comp) != 0
		return result
	else
		return kMultvec(kb(id, 0.0))
	end
end
#########################################################
function Base.:+(A::kMultvec,B::kMultvec)
    return mvsum(A,B)
end
function Base.:-(A::kMultvec,B::kMultvec)
	l=length(B.comp)
	if l != 0
		for i=1:l
			B.comp[i].scl = -B.comp[i].scl
		end
	else
		B = kMultvec([kb(id, false, 0.0)])
	end
	return mvsum(A,B)
end
#########################################################
"""
```
reduct(A::kMultvec)
reduct(A::pMultvec)
reduct(A::cMultvec)

```
Auxiliary function.
"""
function reduct(A::kMultvec)
	l = length(A.comp)
	if l != 0
		a = kMultvec([A.comp[1]])
		B = Vector{kb}(l-1)
		for i=2:l
			B[i-1] = A.comp[i]
		end
		C = mvsum(a, kMultvec(B))
		return C
	else
		return kMultvec([kb(id, 0.0)])
	end
end
#########################################################
function geoprod(A::kMultvec, B::kMultvec)
	l = length(A.comp)
	m = length(B.comp)
	AB = Matrix{kb}(l,m)
	for i=1:l
		for j=1:m
			AB[i,j] = geoprod(A.comp[i], B.comp[j])
		end
	end
	p = kMultvec(reshape(AB,l*m))
	p = reduct(p)
	result = remove(p, 0.0)
	if length(result.comp) != 0
		return result
	else
		return kMultvec([kb(id, 0.0)])
	end
end
#########################################################
function Base.:∘(A::kMultvec,B::kMultvec)
    return geoprod(A,B)
end
#########################################################
function inner(A::kMultvec, B::kMultvec)
	l = length(A.comp)
	m = length(B.comp)
	AB = Matrix{kb}(l,m)
	for i=1:l
		for j=1:m
			AB[i,j] = inner(A.comp[i], B.comp[j])
		end
	end
	p = kMultvec(reshape(AB,l*m))
	p = reduct(p)
	result = remove(p, 0)
	if length(result.comp) != 0
		return result
	else
		return kMultvec([kb(id, 0.0)])
	end
end
#########################################################
function Base.:⋅(A::kMultvec,B::kMultvec)
    return inner(A,B)
end
#########################################################
function outer(A::kMultvec, B::kMultvec)
	l = length(A.comp)
	m = length(B.comp)
	AB = Matrix{kb}(l,m)
	for i=1:l
		for j=1:m
			AB[i,j] = outer(A.comp[i], B.comp[j])
		end
	end
	p = kMultvec(reshape(AB,l*m))
	p = reduct(p)
	result = remove(p, 0)
	if length(result.comp) != 0
		return result
	else
		return kMultvec([kb(id, 0.0)])
	end
end
#########################################################
function Base.:^(A::kMultvec,B::kMultvec)
    return outer(A,B)
end
#########################################################
function copy(A::kMultvec)
	l = length(A.comp)
	X = Vector{kb}(l)
	for i=1:l
		X[i] = copy(A.comp[i])
	end
	return kMultvec(X)
end
#########################################################
function dual(A::kMultvec)
	l = length(A.comp)
	d = Vector{kb}(l)
	for i=1:l
		d[i] = dual(A.comp[i])
	end
	return kMultvec(d)
end
#########################################################
function scalar(A::kMultvec, B::kMultvec)
	l = length(A.comp)
	m = length(B.comp)
	AB = Matrix{kb}(l,m)
	for i=1:l
		for j=1:m
			AB[i,j] = bscalar(A.comp[i], B.comp[j])
		end
	end
	p = kMultvec(reshape(AB,l*m))
	p = reduct(p)
	result = remove(p, 0)
	if length(result.comp) != 0
		return result.comp[1].scl
	else
		return 0.0
	end
end
#########################################################
function Base.:*(A::kMultvec,B::kMultvec)
    return scalar(A,B)
end
#########################################################
function mvreverse(A::kMultvec)
	l = length(A.comp)
	B = copy(A)
	for i=1:l
		B.comp[i] = mvreverse(B.comp[i])
	end
	return B
end
#########################################################
function magnitude(A::kMultvec)
	B = mvreverse(A)
	return sqrt(scalar(A, B))
end
#########################################################
"""
```
bltomv(A::kblade)
bltomv(A::pblade)
bltomv(A::cblade)

```
Auxiliary function.
"""
function bltomv(A::kblade)
	l = length(A.conj)
	X = A.conj[1]
	for i=2:l
		X = outer(X, A.conj[i])
	end
	return X
end
#########################################################
function geoprod(A::kblade, B::kblade)
	X = bltomv(A)
	Y = bltomv(B)
	result = geoprod(X, Y)
	return result
end
#########################################################
function Base.:∘(A::kblade,B::kblade)
    return geoprod(A,B)
end
#########################################################
function outer(A::kblade, B::kblade)
	X = bltomv(A)
	Y = bltomv(B)
	result = outer(X, Y)
	return result
end
#########################################################
function Base.:^(A::kblade,B::kblade)
    return outer(A,B)
end
#########################################################
function inner(A::kblade, B::kblade)
	X = bltomv(A)
	Y = bltomv(B)
	result = inner(X, Y)
	return result
end
#########################################################
function Base.:⋅(A::kblade,B::kblade)
    return inner(A,B)
end
#########################################################
function scalar(A::kblade, B::kblade)
	X = bltomv(A)
	Y = bltomv(B)
	result = scalar(X, Y)
	return result
end
#########################################################
function Base.:*(A::kblade,B::kblade)
    return scalar(A,B)
end
#########################################################
function mvsum(A::kblade, B::kblade)
	AA = bltomv(A)
	BB = bltomv(B)
	return mvsum(AA,BB)
end
#########################################################
function Base.:+(A::kblade,B::kblade)
    return mvsum(A,B)
end
function Base.:-(A::kblade,B::kblade)
	AA = bltomv(A)
	BB = bltomv(B)
	l = length(BB.comp)
	for i=1:l
		BB.comp[i].scl = -BB.comp[i].scl
	end
    return mvsum(AA,BB)
end
#########################################################
function mvreverse(A::kblade)
	l = length(A.conj)
	T = Vector{kMultvec}(l)
	fill!(T, kMultvec(kb(id, 0.0)))
	for i=1:l
		T[i] = A.conj[l-i+1]
	end
	return kblade(T)
end
#########################################################
function magnitude(A::kblade)
	X = scalar(A, mvreverse(A))
	return sqrt(X)
end
#########################################################
function dual(A::kblade)
	X = bltomv(A)
	return dual(X)
end
#########################################################
function copy(A::kblade)
	l = length(A.conj)
	X = Vector{kMultvec}(l)
	for i=1:l
		X[i] = copy(A.conj[i])
	end
	return kblade(X)
end
#########################################################
"""
```
inverse(A::kblade)
inverse(A::pblade)
inverse(A::cblade)

```
Returns the inverse of A if A is a non-null-blade or the pseudoinverse A if A is a null-blade.
"""
function inverse(A::kblade)
	div = scalar(A, mvreverse(A))
	l = length(A.conj)
	aux = mvreverse(A)
	X = copy(aux)
	m = length(aux.conj[1].comp)
	for j=1:m
		X.conj[1].comp[j].scl = (aux.conj[1].comp[j].scl)*(div)^(-1)
	end
	return X
end
#########################################################
"""
```
projection(A::kblade, N::kblade)
projection(A::pblade, N::pblade)

```
Returns the projection of the blade A onto the blade N.
"""
function projection(A::kblade, N::kblade)
	N2 = inverse(N)
	result = geoprod(inner(A, N2), bltomv(N))
	return result
end
#########################################################
"""
```
rejection(A::kblade, N::kblade)
rejection(A::pblade, N::pblade)

```
Returns the rejection of the blade A from the blade N.
"""
function rejection(A::kblade, N::kblade)
	aux = projection(A, N)
	l = length(aux.comp)
	for i=1:l
		aux.comp[i].scl = -(aux.comp[i].scl)
	end
	result = mvsum(bltomv(A), aux)
	return result
end
#########################################################
function show(io::IO, a::kb)
	if (a.scl == -1) && (grade(a) == 0)
		print(io, "-1")
	elseif (a.scl == -1) && (grade(a) != 0)
		print(io, "-")
	elseif (a.scl != 1) || (grade(a) == 0)
		print(io, "$(a.scl)")
	end
	flag = 0
	for i=1:length(a.e)
		if a.e[i] == true
			flag = 1
		end
	end
	if flag == 1 && a.scl != 0
		print(io, "e")
	end
	for i=1:length(a.e)
		if a.e[i] == true && a.scl != 0
			print(io, "$i")
		end
	end
end
#########################################################
function show(io::IO, A::kMultvec)
	l = length(A.comp)
	if l != 0
		for i=1:(l-1)
			print(io, "$(A.comp[i]) + ")
		end
		print(io, "$(A.comp[l])")
	else
		print(io, 0)
	end
end
#########################################################
function show(io::IO, X::kblade)
	l = length(X.conj)
	for i=1:l-1
		print(io, "($(X.conj[i]))∧")
	end
	print(io, "($(X.conj[l]))")
end
#########################################################
"""
```
pb(a::Vector{Bool}, b::Bool, c::Number)
```
Creates a scaled basis element of geometric algebra space from a, b and c.
The b::Bool parameter represent the conponet of the basis that squares to -1.

"""
type pb
	eb::Vector{Bool}
	ep::Bool
	scl::Number

	pb(a::Vector{Bool}, b::Bool, c::Number) = new(a, b, c)
	pb(a::Vector{Bool}, b::Bool) = new(a, b, 1.0)
	pb(b::Number) = new(id, false, b)
end
#########################################################
"""
```
pMultvec(X::Vector{pb})
```

Creates a multivector in the geometric algebra (Gp,1) composed by the sum of the coordinates of X.

"""
type pMultvec
	comp::Vector{pb}

	pMultvec(a::pb) = new([a])
	pMultvec(X::Vector{pb}) = new(X)
end
#########################################################
function grade(a::pb)
	l = length(a.eb)
	cont = 0
	for i=1:l
		if a.eb[i] == true
			cont = cont + 1
		end
	end
	if a.ep == true
		cont = cont + 1
	end
	return cont
end
#########################################################
function show(io::IO, a::pb)
	l = length(a.eb)
	print(io, "$(a.scl)")
	if grade(a) != 0 && a.scl != 0
		print(io, "e")
		for i=1:l-1
			if a.eb[i] == true
				print(io, "$i")
			end
		end
		if a.eb[l] == true
			print(io, "₊")
		end
		if a.ep == true
			print(io, "₋")
		end
	end
end
#########################################################
function show(io::IO, A::pMultvec)
	l = length(A.comp)
	if l == 0
		print(io, "0")
	else
		print(io, "$(A.comp[1])")
		for i=2:l
			if A.comp[i].scl != 0
				print(io, " + $(A.comp[i])")
			else
				print(io, " + 0.0")
			end
		end
	end
end
#########################################################
"""
```
kbtopb(a::kb)
```
Auxiliary function.
"""
function kbtopb(a::kb)
	b = pb(a.e, false, a.scl)
	return b
end
#########################################################
function copy(a::pb)
	eb = copy(a.eb)
	ep = copy(a.ep)
	esc = copy(a.scl)
	result = pb(eb,ep, esc)
	return result
end
#########################################################
function geoprod(a::pb, b::pb)
	l = length(a.eb)
	if a.ep == b.ep == true
		aux = -1
		ep = false
	elseif a.ep == b.ep == false
		aux = 1
		ep = false
	else
		aux = 1
		ep = true
	end
	c = Vector{Bool}(l+1)
	d = Vector{Bool}(l+1)
	for i=1:l
		c[i] = a.eb[i]
		d[i] = b.eb[i]
	end
	c[l+1] = a.ep
	d[l+1] = b.ep
	f = geoprod(kb(c, a.scl), kb(d, b.scl))
	scl = aux*f.scl
	g = Vector{Bool}(l)
	for i=1:l
		g[i] = f.e[i]
	end
	result = pb(g, ep, scl)
	return result
end
#########################################################
function Base.:∘(a::pb,b::pb)
    return geoprod(A,B)
end
#########################################################
function inner(a::pb, b::pb)
	if grade(geoprod(a,b)) == abs(grade(a) - grade(b))
		return geoprod(a,b)
	else
		x = Vector{Bool}(length(a.eb))
		fill!(x, false)
		return pb(x, false, 0.0)
	end
end
#########################################################
function Base.:⋅(a::pb,b::pb)
    return inner(A,B)
end
#########################################################
function outer(a::pb, b::pb)
	if grade(geoprod(a,b)) == (grade(a) + grade(b))
		return geoprod(a,b)
	else
		x = Vector{Bool}(length(a.eb))
		fill!(x, false)
		return pb(x, false, 0.0)
	end
end
#########################################################
function Base.:^(a::pb,b::pb)
    return outer(A,B)
end
#########################################################
function scalar(a::pb, b::pb)
	if grade(geoprod(a,b)) == 0
		return geoprod(a,b).scl
	else
		x = Vector{Bool}(length(a.eb))
		fill!(x, false)
		return pb(x, false, 0.0)
	end
end
#########################################################
function bscalar(a::pb, b::pb)
	if grade(geoprod(a,b)) == 0
		return geoprod(a,b)
	else
		x = Vector{Bool}(length(a.eb))
		fill!(x, false)
		return pb(x, false, 0.0)
	end
end
#########################################################
function Base.:*(a::pb,b::pb)
    return scalar(A,B)
end
#########################################################
function mvsum(a::pb, b::pb)
	l = length(a.eb)
	eb = Vector{Bool}(l)
	fill!(eb, false)
	if a.eb == b.eb && a.ep == b.ep
		scl = a.scl + b.scl
	 	return pb(a.eb, a.ep, scl)
	else
		return pMultvec([a,b])
	end
end
#########################################################
function Base.:+(a::pb,b::pb)
    return mvsum(a,b)
end
function Base.:-(a::pb,b::pb)
	b.scl = -b.scl
	return mvsum(a,b)
end
#########################################################
function mvreverse(a::pb)
	k = grade(a)
	eb = copy(a.eb)
	ep = copy(a.ep)
	b = pb(eb,ep, 1)
	b.scl = (-1)^(k*(k-1)/2)*a.scl
	return b
end
#########################################################
"""
```
conjugate(a::pb)
conjugate(X::pMultvec)
conjugate(A::pblade)
conjugate(A::cblade)
```
Returns the basis element a, multivector X or blade A.
"""
function conjugate(a::pb)
	grm = 0
	if a.ep == true
		grm = 1
	end
	b = mvreverse(a)
	b.scl = b.scl*(-1)^grm
	return b
end
#########################################################
function dual(a::pb)
	l = length(a.eb)
	eb = Vector{Bool}(l)
	fill!(eb, true)
	I = pb(eb,true,1.0)
	Inv = mvreverse(I)
	dual = geoprod(a, Inv)
	return dual
end
#########################################################
"""
```
prtore(a::pb)
```
Auxiliary function.
"""
function prtore(a::pb)
	l = length(a.eb)
	e = Vector{Bool}(l+1)
	for i=1:l
		e[i] = a.eb[i]
	end
	e[l+1] = a.ep
	b = kb(e, a.scl)
	return b
end
#########################################################
"""
```
retopr(a::kb)
```
Auxiliary function.
"""
function retopr(a::kb)
	l = length(a.e)
	eb = Vector{Bool}(l-1)
	ep = a.e[l]
	for i=1:l-1
		eb[i] = a.e[i]
	end
	b = pb(eb,ep,a.scl)
end
#########################################################
"""
```
plen(A::pMultvec)
plen(A::pblade)
```
Auxiliary function.
"""
function plen(A::pMultvec)
	l = length(A.comp)
	if l != 0
		r = length(A.comp[1].eb)
		eb = fill!(Vector{Bool}(r), false)
		return eb
	else
		return id
	end
end
#########################################################
"""
```
mvtopmv(A::kMultvec)
```
Auxiliary function.
"""
function mvtopmv(A::kMultvec)
	l = length(A.comp)
	X = Vector{pb}(l)
	for i=1:l
		X[i] = kbtopb(A.comp[i])
	end
	return pMultvec(X)
end
#########################################################
function remove(A::pMultvec, b::Number)
	l = length(A.comp)
	aux = 0
	for i=1:l
		if A.comp[i].scl == b
			aux = aux+1
		end
	end
	B = Vector{pb}(l-aux)
	for i=1:l
		for j=1:l-1
			if A.comp[j].scl == b
				aux2 = A.comp[j+1]
				A.comp[j+1] = A.comp[j]
				A.comp[j] = aux2
			end
		end
	end
	for i=1:(l-aux)
		B[i] = A.comp[i]
	end
	return pMultvec(B)
end
#########################################################
function copy(A::pMultvec)
	l = length(A.comp)
	X = Vector{pb}(l)
	for i=1:l
		X[i] = copy(A.comp[i])
	end
	return pMultvec(X)
end
#########################################################
function mvsum(A::pMultvec, B::pMultvec)
	l = length(A.comp)
	m = length(B.comp)
	A2 = Vector{kb}(l)
	B2 = Vector{kb}(m)
	for i=1:l
		A2[i] = prtore(A.comp[i])
	end
	for i=1:m
		B2[i] = prtore(B.comp[i])
	end
	A2 = kMultvec(A2)
	B2 = kMultvec(B2)
	AB2 = mvsum(A2,B2)
	n = length(AB2.comp)
	AB = Vector{pb}(n)
	for i=1:n
		AB[i] = retopr(AB2.comp[i])
	end
	AB = pMultvec(AB)
	return remove(AB, 0.0)
end
#########################################################
function Base.:+(A::pMultvec,B::pMultvec)
    return mvsum(A,B)
end
function Base.:-(A::pMultvec,B::pMultvec)
	l=length(B.comp)
	if l != 0
		for i=1:l
			B.comp[i].scl = -B.comp[i].scl
		end
	else
		eb = plen(B)
		B = pMultvec([pb(eb, false, 0.0)])
	end
	return mvsum(A,B)
end
#########################################################
function reduct(A::pMultvec)
	l = length(A.comp)
	eb = plen(A)
	if l != 0
		a = pMultvec([A.comp[1]])
		B = Vector{pb}(l-1)
		for i=2:l
			B[i-1] = A.comp[i]
		end
		C = mvsum(a, pMultvec(B))
		return C
	else
		return pMultvec([pb(eb, false, 0.0)])
	end
end
#########################################################
function geoprod(A::pMultvec, B::pMultvec)
	l = length(A.comp)
	m = length(B.comp)
	eb = plen(A)
	AB = Matrix{pb}(l,m)
	for i=1:l
		for j=1:m
			AB[i,j] = geoprod(A.comp[i], B.comp[j])
		end
	end
	p = pMultvec(reshape(AB,l*m))
	p = reduct(p)
	result = remove(p, 0)
	if length(result.comp) != 0
		return result
	else
		return pMultvec([pb(eb, false, 0.0)])
	end
end
#########################################################
function Base.:∘(A::pMultvec,B::pMultvec)
    return geoprod(A,B)
end
#########################################################
function inner(A::pMultvec, B::pMultvec)
	l = length(A.comp)
	m = length(B.comp)
	eb = plen(A)
	AB = Matrix{pb}(l,m)
	for i=1:l
		for j=1:m
			AB[i,j] = inner(A.comp[i], B.comp[j])
		end
	end
	p = pMultvec(reshape(AB,l*m))
	p = reduct(p)
	result = remove(p, 0)
	if length(result.comp) != 0
		return result
	else
		return pMultvec([pb(eb, false, 0.0)])
	end
end
#########################################################
function Base.:⋅(A::pMultvec,B::pMultvec)
    return inner(A,B)
end
#########################################################
function outer(A::pMultvec, B::pMultvec)
	l = length(A.comp)
	m = length(B.comp)
	AB = Matrix{pb}(l,m)
	eb = plen(A)
	for i=1:l
		for j=1:m
			AB[i,j] = outer(A.comp[i], B.comp[j])
		end
	end
	p = pMultvec(reshape(AB,l*m))
	p = reduct(p)
	result = remove(p, 0)
	if length(result.comp) != 0
		return result
	else
		return pMultvec([pb(eb, false, 0.0)])
	end
end
#########################################################
function Base.:^(A::pMultvec,B::pMultvec)
    return outer(A,B)
end
#########################################################
function scalar(A::pMultvec, B::pMultvec)
	l = length(A.comp)
	m = length(B.comp)
	AB = Matrix{pb}(l,m)
	for i=1:l
		for j=1:m
			AB[i,j] = bscalar(A.comp[i], B.comp[j])
		end
	end
	p = pMultvec(reshape(AB,l*m))
	p = reduct(p)
	result = remove(p, 0)
	if length(result.comp) != 0
		return result.comp[1].scl
	else
		return 0.0
	end
end
#########################################################
function Base.:*(A::pMultvec,B::pMultvec)
    return scalar(A,B)
end
#########################################################
function conjugate(A::pMultvec)
	l = length(A.comp)
	B = Vector{pb}(l)
	for i=1:l
		B[i] = conjugate(A.comp[i])
	end
	return pMultvec(B)
end
#########################################################
function magnitude(A::pMultvec)
	scl = scalar(A,conjugate(A))
	result = sqrt(scl)
	return result
end
#########################################################
function dual(A::pMultvec)
	k = length(A.comp)
	d = Vector{pb}(k)
	for i=1:k
		d[i] = dual(A.comp[i])
	end
	return pMultvec(d)
end
#########################################################
"""
```
mvectovec(X::Vector{pMultvec})
```
Auxiliary function.
"""
function mvectovec(X::Vector{pMultvec})
	l = length(X)
	flag = 0
	for i=1:l
		m = length(X[i].comp)
		for j=1:m
			if grade(X[i].comp[j]) != 1
				flag = 1
			end
		end
	end
	if flag == 1
		error("Grade error (mvectovec function)")
		return 0
	else
		n = length(X[1].comp[1].eb)
		A = Matrix{Number}(l, n+1)
		fill!(A, 0.0)
		for i=1:l
			m = length(X[i].comp)
			for j=1:m
				for k=1:n
					if X[i].comp[j].eb[k] == true
						A[i,k] = X[i].comp[j].scl
					end
				end
				if X[i].comp[j].ep == true
					A[i,n+1] = X[i].comp[j].scl
				end
			end
		end
		return A
	end
end
#########################################################
"""
```
pblade(X::Vector{pMultvec})
```
Creates a blade in the geometric algebra (Gp,1) composed by the outer product of the elements of X.

The blade is created only if the multivectors that compound X are both 1-vectors and X is a L.i set

"""
type pblade
	conj::Vector{pMultvec}
	function pblade(X::Vector{pMultvec})
		T = mvectovec(X)
		if det(T*transpose(T)) != 0
			new(X)
		else
			error("L.d(pblade conversion)")
		end
	end
end
#########################################################
function show(io::IO, A::pblade)
	if length(A.conj) != 0
		print(io, "($(A.conj[1]))")
		for i=2:length(A.conj)
			print(io, "∧($(A.conj[i]))")
		end
	else
		print(io, "0")
	end
end
#########################################################
function bltomv(A::pblade)
	l = length(A.conj)
	X = A.conj[1]
	for i=2:l
		X = outer(X, A.conj[i])
	end
	return X
end
#########################################################
function geoprod(A::pblade, B::pblade)
	X = bltomv(A)
	Y = bltomv(B)
	result = geoprod(X, Y)
	return result
end
#########################################################
function Base.:∘(A::pblade, B::pblade)
    return geoprod(A,B)
end
#########################################################
function outer(A::pblade, B::pblade)
	X = bltomv(A)
	Y = bltomv(B)
	result = outer(X, Y)
	return result
end
#########################################################
function Base.:^(A::pblade, B::pblade)
    return outer(A,B)
end
#########################################################
function inner(A::pblade, B::pblade)
	X = bltomv(A)
	Y = bltomv(B)
	result = inner(X, Y)
	return result
end
#########################################################
function Base.:⋅(A::pblade, B::pblade)
    return inner(A,B)
end
#########################################################
function scalar(A::pblade, B::pblade)
	X = bltomv(A)
	Y = bltomv(B)
	result = scalar(X, Y)
	return result
end
#########################################################
function Base.:*(A::pblade, B::pblade)
    return scalar(A,B)
end
#########################################################
function plen(A::pblade)
	B = bltomv(A)
	l = length(B.comp)
	if l != 0
		r = length(B.comp[1].eb)
		eb = fill!(Vector{Bool}(r), false)
		return eb
	else
		return id
	end
end
#########################################################
function mvreverse(A::pblade)
	l = length(A.conj)
	eb = plen(A)
	T = Vector{pMultvec}(l)
	fill!(T, pMultvec([pb(eb, false, 0.0)]))
	for i=1:l
		T[i] = A.conj[l-i+1]
	end
	return pblade(T)
end
#########################################################
function conjugate(A::pblade)
	l = length(A.conj)
	T = Vector{pMultvec}(l)
	eb = plen(A)
	fill!(T, pMultvec([pb(eb, false, 0.0)]))
	for i=1:l
		T[i] = conjugate(A.conj[l-i+1])
	end
	return pblade(T)
end
#########################################################
function magnitude(A::pblade)
	X = scalar(A, conjugate(A))
	return sqrt(X)
end
#########################################################
function dual(A::pblade)
	X = bltomv(A)
	return dual(X)
end
#########################################################
function copy(A::pblade)
	l = length(A.conj)
	X = Vector{pMultvec}(l)
	for i=1:l
		X[i] = copy(A.conj[i])
	end
	return pblade(X)
end
#########################################################
function inverse(A::pblade)
	if geoprod(A,A).comp[1].scl != 0
		div = scalar(A, mvreverse(A))
		l = length(A.conj)
		aux = mvreverse(A)
		X = copy(aux)
		m = length(aux.conj[1].comp)
		for j=1:m
			X.conj[1].comp[j].scl = (aux.conj[1].comp[j].scl)*(div)^(-1)
		end
	else
		div = inner(A, conjugate(A)).comp[1].scl
		l = length(A.conj)
		aux = conjugate(A)
		X = copy(aux)
		m = length(aux.conj[1].comp)
		for j=1:m
			X.conj[1].comp[j].scl = (aux.conj[1].comp[j].scl)*(div)^(-1)
		end
	end
	return X
end
#########################################################
function projection(A::pblade, N::pblade)
	X = inner(inner(A, inverse(N)), bltomv(N))
	return X
end
#########################################################
function rejection(A::pblade, N::pblade)
	X = projection(A, N)
	l = length(X.comp)
	for i=1:l
		X.comp[i].scl = (-1)*(X.comp[i].scl)
	end
	result = mvsum(bltomv(A), X)
	return result
end
#########################################################
"""
```
cb(a::Vector{Bool},b::Bool,c::Bool,d::Number)
```
Creates a scaled basis element of geometric algebra space (Gp+1,1) from a, b, c and d using e∞ and e∘ as basis elements instead of the regular basis.

"""
type cb
	er::Vector{Bool}
	ei::Bool
	eo::Bool
	scl::Number
	function cb(a::Vector{Bool},b::Bool,c::Bool,d::Number)
		if b == c == false
			new(a,b,c,d)
		elseif b != c
			new(a,b,c,d)
		elseif b == c
			new(a, false, false, -1.0*d)
		end
	end
end
#########################################################
"""
```
cMultvec(X::Vector{cb})
```

Creates a multivector in the geometric algebra (Gp+1,1) composed by the sum of the coordinates of X.

"""
type cMultvec
	comp::Vector{cb}

end
#########################################################
function grade(a::cb)
	l = length(a.er)
	cont = 0
	for i=1:l
		if a.er[i] == true
			cont = cont + 1
		end
	end
	if a.ei == true
		cont = cont + 1
	end
	if a.eo == true
		cont = cont + 1
	end
	return cont
end
#########################################################
function show(io::IO, a::cb)
	l = length(a.er)
	print(io, "$(a.scl)")
	if grade(a) != 0
		print(io, "e")
		for i=1:l
			if a.er[i] == true
				print(io, "$i")
			end
		end
		if a.ei == true
			print(io, "∞")
		end
		if a.eo == true
			print(io, "ₒ")
		end
	end
end
#########################################################
function show(io::IO, A::cMultvec)
	l = length(A.comp)
	if l == 0
		print(io, "0")
	else
		print(io, "$(A.comp[1])")
		for i=2:l
			if A.comp[i].scl != 0
				print(io, " + $(A.comp[i])")
			else
				print(io, " + 0.0")
			end
		end
	end
end
#########################################################
function copy(A::cMultvec)
	l = length(A.comp)
	X = Vector{cb}(l)
	for i=1:l
		X[i] = copy(A.comp[i])
	end
	return cMultvec(X)
end
#########################################################
function geoprod(a::cb, b::cb)
	l=length(a.er)
	aux = 1
	a1 = kb(a.er, a.scl)
	b1 = kb(b.er, b.scl)
	ab1 = geoprod(a1, b1)
	ei = false
	eo = false
	if a.ei == b.ei == true || a.eo == b.eo == true
		aux = 0
		fill!(ab1.e,false)
	elseif a.ei == b.eo == true || a.eo == b.ei == true
		aux = -1
		ei = false
		eo = false
	end
	if a.ei == true && a.eo == b.ei == b.eo == false || b.ei == true && a.ei == a.eo == b.eo == false
		ei = true
		eo = false
	end
	if a.eo == true && a.ei == b.ei == b.eo == false || b.eo == true && a.ei == a.eo == b.ei == false
		ei = false
		eo = true
	end
	ab = cb(ab1.e, ei, eo, aux*ab1.scl)
	return ab
end
#########################################################
function Base.:∘(a::cb, b::cb)
    return geoprod(a,b)
end
#########################################################
function inner(a::cb, b::cb)
	if grade(geoprod(a,b)) == abs(grade(a) - grade(b))
		return geoprod(a,b)
	else
		return cb(id, false, false, 0.0)
	end
end
#########################################################
function Base.:⋅(a::cb, b::cb)
    return inner(a,b)
end
#########################################################
function outer(a::cb, b::cb)
	if grade(geoprod(a,b)) == grade(a) + grade(b)
		return geoprod(a,b)
	else
		return cb(id, false, false, 0.0)
	end
end
#########################################################
function Base.:^(a::cb, b::cb)
    return outer(a,b)
end
#########################################################
function scalar(a::cb, b::cb)
	if grade(geoprod(a,b)) == 0
		return geoprod(a,b).scl
	else
		return cb(id, false, false, 0.0)
	end
end
#########################################################
function bscalar(a::cb, b::cb)
	if grade(geoprod(a,b)) == 0
		return geoprod(a,b)
	else
		return cb(id, false, false, 0.0)
	end
end
#########################################################
function Base.:*(a::cb, b::cb)
    return scalar(a,b)
end
#########################################################
function mvsum(a::cb, b::cb)
	l = length(a.er)
	er = Vector{Bool}(l)
	fill!(er, false)
	if a.er == b.er && a.ei == b.ei && a.eo == b.eo
		scl = a.scl + b.scl
	 	return cb(a.er, a.ei,a.eo, scl)
	else
		return cMultvec([a,b])
	end
end
#########################################################
function Base.:+(a::cb, b::cb)
    return mvsum(a,b)
end
function Base.:-(a::cb, b::cb)
	b.scl = - b.scl
    return mvsum(a,b)
end
#########################################################
function mvreverse(a::cb)
	k = grade(a)
	er = copy(a.er)
	ei = copy(a.ei)
	eo = copy(a.eo)
	b = cb(er, ei, eo, 1.0)
	b.scl = (-1)^(k*(k-1)/2)*a.scl
	return b
end
#########################################################
"""
```
cbtopb(a::cb)
```
Auxiliary function.
"""
function cbtopb(a::cb)
	l = length(a.er)
	ea = fill!(Vector{Bool}(l+1), false)
	eb = copy(ea)
	eb[l+1] = true
	X = pMultvec(pb(ea, false, 0.0))
	Y = copy(X)
	f = 0
	er = Vector{Bool}(l+1)
	for i=1:l
		er[i] = a.er[i]
	end
	er[l+1] = false
	if a.ei == true
		X = pMultvec([pb(ea, true, 1.0), pb(eb, false, 1.0)])
		f = 1
	end
	if a.eo == true
		Y = pMultvec([pb(ea, true, 0.5), pb(eb, false, -0.5)])
		f = 1
	end
	if f == 1
		p1 = geoprod(pMultvec(pb(er, false, a.scl)), X)
		p2 = geoprod(pMultvec(pb(er, false, a.scl)), Y)
		return mvsum(p1, p2)
	else
		return pb(er, false, a.scl)
	end
end
#########################################################
"""
```
cbtore(a::cb)
```
Auxiliary function.
"""
function cbtore(a::cb)
	l = length(a.er)
	e = Vector{Bool}(l+2)
	for i=1:l
		e[i] = a.er[i]
	end
	e[l+1] = a.ei
	e[l+2] = a.eo
	b = kb(e, a.scl)
	return b
end
#########################################################
"""
```
retocb(a::kb)
```
Auxiliary function.
"""
function retocb(a::kb)
	l = length(a.e)
	er = Vector{Bool}(l-2)
	for i=1:l-2
		er[i] = a.e[i]
	end
	ei = a.e[l-1]
	eo = a.e[l]
	b = cb(er, ei, eo, a.scl)
	return b
end
#########################################################
function remove(A::cMultvec, b::Number)
	l = length(A.comp)
	aux = 0
	for i=1:l
		if A.comp[i].scl == b
			aux = aux+1
		end
	end
	B = Vector{cb}(l-aux)
	if aux != 0
		for i=1:l
			for j=1:l-1
				if A.comp[j].scl == b
					aux2 = A.comp[j+1]
					A.comp[j+1] = A.comp[j]
					A.comp[j] = aux2
				end
			end
		end
		for i=1:(l-aux)
			B[i] = A.comp[i]
		end
		return cMultvec(B)
	else
		return A
	end
end
#########################################################
function mvsum(A::cMultvec, B::cMultvec)
	l = length(A.comp)
	m = length(B.comp)
	A2 = Vector{kb}(l)
	B2 = Vector{kb}(m)
	for i=1:l
		A2[i] = cbtore(A.comp[i])
	end
	for i=1:m
		B2[i] = cbtore(B.comp[i])
	end
	A2 = kMultvec(A2)
	B2 = kMultvec(B2)
	AB2 = mvsum(A2,B2)
	n = length(AB2.comp)
	AB = Vector{cb}(n)
	for i=1:n
		AB[i] = retocb(AB2.comp[i])
	end
	AB = cMultvec(AB)
	return remove(AB, 0.0)
end
#########################################################
function Base.:+(A::cMultvec,B::cMultvec)
    return mvsum(A,B)
end
function Base.:-(A::cMultvec,B::cMultvec)
	l=length(B.comp)
	if l != 0
		for i=1:l
			B.comp[i].scl = -B.comp[i].scl
		end
	else
		B = cMultvec([cb(id,false, false, 0.0)])
	end
	return mvsum(A,B)
end
#########################################################
function reduct(A::cMultvec)
	l = length(A.comp)
	if l != 0
		a = cMultvec([A.comp[1]])
		B = Vector{cb}(l-1)
		for i=2:l
			B[i-1] = A.comp[i]
		end
		C = mvsum(a, cMultvec(B))
		return C
	else
		return cMultvec([cb(id,false, false, 0.0)])
	end
end
#########################################################
function geoprod(A::cMultvec, B::cMultvec)
	l = length(A.comp)
	m = length(B.comp)
	AB = Matrix{cb}(l,m)
	for i=1:l
		for j=1:m
			AB[i,j] = geoprod(A.comp[i], B.comp[j])
		end
	end
	p = cMultvec(reshape(AB,l*m))
	p = reduct(p)
	result = remove(p, 0)
	if length(result.comp) != 0
		return result
	else
		return cMultvec([cb(id,false, false, 0.0)])
	end
end
#########################################################
function Base.:∘(A::cMultvec,B::cMultvec)
    return geoprod(A,B)
end
#########################################################
function inner(A::cMultvec, B::cMultvec)
	l = length(A.comp)
	m = length(B.comp)
	AB = Matrix{cb}(l,m)
	for i=1:l
		for j=1:m
			AB[i,j] = inner(A.comp[i], B.comp[j])
		end
	end
	p = cMultvec(reshape(AB,l*m))
	p = reduct(p)
	result = remove(p, 0)
	if length(result.comp) != 0
		return result
	else
		return cMultvec([cb(id,false, false, 0.0)])
	end
end
#########################################################
function Base.:⋅(A::cMultvec,B::cMultvec)
    return inner(A,B)
end
#########################################################
function outer(A::cMultvec, B::cMultvec)
	l = length(A.comp)
	m = length(B.comp)
	AB = Matrix{cb}(l,m)
	for i=1:l
		for j=1:m
			AB[i,j] = outer(A.comp[i], B.comp[j])
		end
	end
	p = cMultvec(reshape(AB,l*m))
	p = reduct(p)
	result = remove(p, 0)
	if length(result.comp) != 0
		return result
	else
		return cMultvec([cb(id,false, false, 0.0)])
	end
end
#########################################################
function Base.:^(A::cMultvec,B::cMultvec)
    return outer(A,B)
end
#########################################################
function scalar(A::cMultvec, B::cMultvec)
	l = length(A.comp)
	m = length(B.comp)
	AB = Matrix{cb}(l,m)
	for i=1:l
		for j=1:m
			AB[i,j] = bscalar(A.comp[i], B.comp[j])
		end
	end
	p = cMultvec(reshape(AB,l*m))
	p = reduct(p)
	result = remove(p, 0)
	if length(result.comp) != 0
		return result.comp[1].scl
	else
		return 0.0
	end
end
#########################################################
function Base.:*(A::cMultvec,B::cMultvec)
    return scalar(A,B)
end
#########################################################
"""
```
pbtocb(a::pb)
```
Auxiliary function.
"""
function pbtocb(a::pb)
	l = length(a.eb)
	eb = fill!(Vector{Bool}(l-1), false)
	for i=1:l-1
		eb[i] = a.eb[i]
	end
	a2 = cb(eb, false, false, a.scl)
	aux1 = cMultvec([cb(id, false, false, 1.0)])
	aux2 = cMultvec([cb(id, false, false, 1.0)])
	f = 0
	if a.eb[l] == true
		aux1 = cMultvec([cb(id, true, false, 0.5), cb(id,false, true, -1.0)])
		f = 1
	end
	if a.ep == true
		aux2 = cMultvec([cb(id, true, false, 0.5), cb(id, false, true, 1.0)])
		f = 1
	end
	b = cb(eb, false, false, a.scl)
	if f == 1
		b = cMultvec([b])
		p1 = geoprod(aux1, aux2)
		return geoprod(b, p1)
	else
		return b
	end
end
#########################################################
function copy(a::cb)
	er = copy(a.er)
	ei = copy(a.ei)
	eo = copy(a.eo)
	scl = copy(a.scl)
	return cb(er,ei,eo,scl)
end
#########################################################
"""
```
pvtocv(A::pMultvec)
```
Auxiliary function.
"""
function pvtocv(A::pMultvec)
	l = length(A.comp)
	B = Vector{Any}(l)
	for i=1:l
		aux = pbtocb(A.comp[i])
		if typeof(aux) == cb
			B[i] = aux
		else
			B[i] = cb(id, false,false, 0.0)
			B = vcat(B, aux.comp)
		end
	end
	B = cMultvec(B)
	B = reduct(B)
	return B
end
#########################################################
"""
```
cvtopv(A::cMultvec)
```
Auxiliary function.
"""
function cvtopv(A::cMultvec)
	l = length(A.comp)
	m = length(A.comp[1].er)
	eb = fill!(Vector{Bool}(m+1), false)
	B = Vector{pb}(l)
	for i=1:l
		aux = cbtopb(A.comp[i])
		if typeof(aux) == pb
			B[i] = aux
		else
			B[i] = pb(eb,false, 0.0)
			B = vcat(B, aux.comp)
		end
	end
	B = pMultvec(B)
	B = reduct(B)
	return pMultvec(B)
end
#########################################################
"""
```
retoaffin(a::Vector{Float64})
```
Takes a vector of Euclidean space into a vector of affine space.
"""
function retoaffin(a::Vector{Float64})
	l = length(a)
	A = Vector{kb}(l+1)
	aux = Vector{Bool}(l+1)
	fill!(aux, false)
	for i=1:l
		aux[i] = true
		A[i] = copy(kb(aux, a[i]))
		aux[i] = false
	end
	aux[l+1] = true
	A[l+1] = copy(kb(aux, 1.0))
	return kMultvec(A)
end
#########################################################
"""
```
affine(A::kMultvec)
```
Takes a vector of Euclidean subspace of Gp into a vector of affine space.
"""
function affine(A::kMultvec)
	X = copy(reduct(A))
	l = length(X.comp)
	m = length(X.comp[1].e)
	f1 = 1
	f2 = 0
	for i=1:l
		if grade(A.comp[i]) != 1
			f1 = 0
		end
		if length(A.comp[i].e[m]) == true
			f2 = 1
		end
	end
	if f1 == f2 == 1
		aux = Vector{Bool}(m)
		fill!(aux, false)
		aux[m] = true
		scl = inner(A, kMultvec(kb(aux, 1.0))).comp[1].scl
		a = copy(X)
		for i=1:l
			X.comp[i].scl = (X.comp[i].scl)/scl
		end
		return X
	elseif f1 == 1 && f2 == 0
		error("Image error")
	elseif f1 == 0 && f2 == 1
		error("Grade error")
	end
end
#########################################################
"""
```
iretoaffin(A::kMultvec)
```
The inverse transformation of "affine" function.

Takes a vector in affine space back into a vector of Euclidean space subspace of Gp.
"""
function iretoaffin(A::kMultvec)
	l = length(A.comp)
	m =	length(A.comp[1].e)
	f1 = 1
	f2 = 0
	aux = 0
	X = reduct(copy(A))
	for i=1:l
		if grade(X.comp[i]) != 1
			f1 = 0
		end
		if length(X.comp[i].e[m]) == true
			f2 = 1
			aux = i
		end
	end
	if f1 == f2 == 1
		X = affine(X)
		e = fill!(Vector{Bool}(m),false)
		X.comp[aux] = kb(e, 0.0)
		return reduct(X)
	elseif f1 == 1 && f2 == 0
		error("Image error")
	elseif f1 == 0 && f2 == 1
		error("Grade error")
	end
end
#########################################################
"""
```
euctoga(x::Vector{Float64})
```
Auxiliary function.
"""
function euctoga(x::Vector{Float64})
	l = length(x)
	X = Vector{kb}(l)
	e = fill!(Vector{Bool}(l), false)
	for i=1:l
		e[i] = true
		X[i] = kb(copy(e), x[i])
		e[i] = false
	end
	return kMultvec(X)
end
#########################################################
"""
```
S(x::Vector{Float64})
```
The stereographic embedding of Euclidean space.

Takes a vector of Euclidean space in its stereographic embedding.
"""
function S(x::Vector{Float64})
	l = length(x)
	x2 = vecdot(x,x)
	aux = 2/(x2 + 1)
	X = Vector{Float64}(l+1)
	for i=1:l
		X[i] = copy(aux*x[i])
	end
	X[l+1] = (x2 - 1)/(x2 + 1)
	return euctoga(X)
end
#########################################################
"""
```
H(x::Vector{Float64})
```
Auxiliary function.
"""
function H(x::Vector{Float64})
	l = length(x)
	y = euctoga(x)
	X = Vector{pb}(l+1)
	for i=1:l
		X[i] = kbtopb(y.comp[i])
	end
	m = length(S.comp[1].e)
	ide = Vector{Bool}(m+1)
	fill!(ide, false)
	X[l+1] = pb(ide, true, 1.0)
	return pMultvec(X)
end
#########################################################
function H(S::kMultvec)
	l = length(S.comp)
	X = Vector{pb}(l+1)
	for i=1:l
		X[i] = kbtopb(S.comp[i])
	end
	m = length(S.comp[1].e)
	ide = Vector{Bool}(m)
	fill!(ide, false)
	X[l+1] = pb(ide, true, 1.0)
	return pMultvec(X)
end
#########################################################
"""
```
pconformal(x::Vector{Float64})
```
Conformal embedding of a Euclidean vector into a vector of Gp+1,1 with new basis elements being e+ and e-.
"""
function pconformal(x::Vector{Float64})
	x2 = vecdot(x,x)
	a = 0.5*(x2 + 1)
	X = H(S(x))
	l = length(X.comp)
	for i=1:l
		X.comp[i].scl = (X.comp[i].scl)*a
	end
	return X
end
function pconformal(x::kMultvec)
	y = mvectovec(x)
	return pconformal(x)
end
#########################################################
"""
```
ipconformal(X::pMultvec)
```
The inverse of "pconformal" function.
"""
function ipconformal(X::pMultvec)
	l=length(X.comp)
	Y = copy(X)
	if l != 0
		m = length(X.comp[1].eb)
	end
	eb = Vector{Bool}(m)
	fill!(eb,false)
	eb2 = copy(eb)
	eb2[m] = true
	ei = pMultvec([pb(eb, true, 1.0), pb(eb2, false, 1.0)])
	eo = pMultvec([pb(eb, true, 0.5), pb(eb2, false, -0.5)])
	d = inner(X, ei)
	n = length(d.comp)
	if n != 0
		div = -d.comp[1].scl
	end
	for i=1:l
		Y.comp[i].scl = Y.comp[i].scl/(div)
	end
	B = pblade([Y])
	result = rejection(B, pblade([ei, eo]))
	return result
end
#########################################################
"""
```
conformal(x::Vector{Float64})
```
Conformal embedding of a Euclidean vector into a vector of Gp+1,1 with new basis elements being e∞ and e∘.
"""
function conformal(x::Vector{Float64})
	x2 = vecdot(x,x)
	a = 0.5*(x2 + 1)
	X = H(S(x))
	l = length(X.comp)
	for i=1:l
		X.comp[i].scl = (X.comp[i].scl)*a
	end
	return pvtocv(X)
end
#########################################################
"""
```
iconformal(X::cMultvec)
```
The inverse of "conformal" function.
"""
function iconformal(X::cMultvec)
	Y = cvtopv(X)
	Y = ipconformal(Y)
	return pvtocv(Y)
end
#########################################################
"""
```
mvectovec(X::Vector{cMultvec})
```
Auxiliary function.
"""
function mvectovec(X::Vector{cMultvec})
	l = length(X)
	Y = Vector{pMultvec}(l)
	for i=1:l
		Y[i] = cvtopv(X[i])
	end
	return mvectovec(Y)
end
#########################################################
type cblade
	conj::Vector{cMultvec}
	function cblade(X::Vector{cMultvec})
		T = mvectovec(X)
		if det(T*transpose(T)) != 0
			new(X)
		else
			error("L.d (cblade conversion)")
		end
	end
end
#########################################################
function show(io::IO, A::cblade)
	if length(A.conj) != 0
		print(io, "($(A.conj[1]))")
		for i=2:length(A.conj)
			print(io, "∧($(A.conj[i]))")
		end
	else
		print(io, "0")
	end
end
#########################################################
"""
```
bltomv(A::cblade)
```
Auxiliary function.
"""
function bltomv(A::cblade)
	l = length(A.conj)
	X = A.conj[1]
	for i=2:l
		X = outer(X, A.conj[i])
	end
	return X
end
#########################################################
function geoprod(A::cblade,B::cblade)
	AA = bltomv(A)
	BB = bltomv(B)
	return geoprod(AA,BB)
end
#########################################################
function Base.:∘(A::cblade,B::cblade)
	return geoprod(A,B)
end
#########################################################
function inner(A::cblade,B::cblade)
	AA = bltomv(A)
	BB = bltomv(B)
	return inner(AA,BB)
end
#########################################################
function Base.:⋅(A::cblade,B::cblade)
	return inner(A,B)
end
#########################################################
function outer(A::cblade,B::cblade)
	AA = bltomv(A)
	BB = bltomv(B)
	return outer(AA,BB)
end
#########################################################
function Base.:^(A::cblade,B::cblade)
	return outer(A,B)
end
#########################################################
function scalar(A::cblade,B::cblade)
	AA = bltomv(A)
	BB = bltomv(B)
	return scalar(AA,BB)
end
#########################################################
function Base.:*(A::cblade,B::cblade)
	return scalar(A,B)
end
#########################################################
function mvsum(A::cblade,B::cblade)
	AA = bltomv(A)
	BB = bltomv(B)
	return mvsum(AA,BB)
end
#########################################################
function Base.:+(A::cblade,B::cblade)
	return mvsum(A,B)
end
function Base.:-(A::cblade,B::cblade)
	AA = bltomv(A)
	BB = bltomv(B)
	return AA-BB
end
#########################################################
"""
```
cbltopbl(A::cblade)
```
Auxiliary function.
"""
function cbltopbl(A::cblade)
	l = length(A.conj)
	V = Vector{pMultvec}(l)
	for i=1:l
		V[i] = cvtopv(A.conj[i])
	end
	B = pblade(V)
	return B
end
#########################################################
"""
```
pbltocbl(A::pblade)
```
Auxiliary function.
"""
function pbltocbl(A::pblade)
	l = length(A.conj)
	V = Vector{cMultvec}(l)
	for i=1:l
		V[i] = pvtocv(A.conj[i])
	end
	B = cblade(V)
	return B
end
#########################################################
function mvreverse(A::cblade)
	B = cbltopbl(A)
	B = mvreverse(B)
	return pbltocbl(B)
end
#########################################################
function conjugate(A::cblade)
	B = cbltopbl(A)
	B = conjugate(B)
	return pbltocbl(B)
end
#########################################################
function magnitude(A::cblade)
	X = A * conjugate(A)
	return sqrt(X)
end
#########################################################
function copy(A::cblade)
	l = length(A.conj)
	X = Vector{cMultvec}(l)
	for i=1:l
		X[i] = copy(A.conj[i])
	end
	return cblade(X)
end
#########################################################
function inverse(A::cblade)
	if (A∘A).comp[1].scl != 0
		div = (A * mvreverse(A))
		l = length(A.conj)
		aux = mvreverse(A)
		X = copy(aux)
		m = length(aux.conj[1].comp)
		for j=1:m
			X.conj[1].comp[j].scl = (aux.conj[1].comp[j].scl)*(div)^(-1)
		end
	else
		div = (A ⋅ conjugate(A)).comp[1].scl
		l = length(A.conj)
		aux = conjugate(A)
		X = copy(aux)
		m = length(aux.conj[1].comp)
		for j=1:m
			X.conj[1].comp[j].scl = (aux.conj[1].comp[j].scl)*(div)^(-1)
		end
	end
	return X
end
#########################################################
