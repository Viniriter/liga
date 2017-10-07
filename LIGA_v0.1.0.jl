import Base.show, Base.copy, Base.subtypes
importall Base.Operators
#########################################################
#########################################################
type kbase
    e::Vector{Bool}
	scl::Number

	kbase(a::Vector{Bool}) =   new(a, 1.0)
	kbase(a::Vector{Bool}, b::Number) =   new(a, b)
end
#########################################################
type kMultvec
    comp::Vector{kbase}

	kMultvec(a::kbase) = new([a])
	kMultvec(A::Vector{kbase}) = new(A)
	kMultvec(n::Int) = new(Vector{kbase}(n))
end
#########################################################
function grade(x::kbase)
	acu = 0
	for i=1:length(x.e)
		if x.e[i] == true
			acu = acu + 1
		end
	end
	return acu
end
#########################################################
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
function copy(a::kbase)
	e = copy(a.e)
	esc = copy(a.scl)
	result = kbase(e, esc)
	return result
end
#########################################################
function bgeoprod(a::kbase, b::kbase)
	l = length(a.e)
	Ea = copy(a.e)
	Eb = copy(b.e)
	ab = kbase(copy(Ea), 0)
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
function Base.:∘(a::kbase,b::kbase)
    return bgeoprod(a,b)
end
#########################################################
function innprod(a::kbase, b::kbase)
	if grade(bgeoprod(a,b)) == abs(grade(a) - grade(b))
		return bgeoprod(a,b)
	else
		x = Vector{Bool}(length(a.e))
		fill!(x, false)
		return kbase(x, 0)
	end
end
#########################################################
function Base.:⋅(a::kbase,b::kbase)
    return innprod(a,b)
end
#########################################################
function outprod(a::kbase, b::kbase)
	if grade(bgeoprod(a,b)) == (grade(a) + grade(b))
		return bgeoprod(a,b)
	else
		x = Vector{Bool}(length(a.e))
		fill!(x, false)
		return kbase(x, 0)
	end
end
#########################################################
function Base.:^(a::kbase,b::kbase)
    return outprod(a,b)
end
#########################################################
function escalarb(a::kbase, b::kbase)
	if grade(bgeoprod(a,b)) == 0
		return bgeoprod(a,b)
	else
		x = Vector{Bool}(length(a.e))
		fill!(x, false)
		return kbase(x, 0)
	end
end
#########################################################
function Base.:*(a::kbase,b::kbase)
    return escalarb(a,b)
end
#########################################################
function bsum(a::kbase, b::kbase)
	l = length(a.e)
	e = Vector{Bool}(l)
	fill!(e, false)
	sume = kbase(e, 0)
	if a.e == b.e
		sume.e = a.e
		sume.scl = a.scl + b.scl
	elseif a.e != b.e
		sume = kMultvec([a,b])
	end
	return sume
end
#########################################################
function Base.:+(a::kbase,b::kbase)
    return bsum(a,b)
end
function Base.:-(a::kbase,b::kbase)
	b.scl = -b.scl
	return bsum(a,b)
end
#########################################################
function breverse(a::kbase)
	k = grade(a)
	e = copy(a.e)
	b = kbase(e)
	b.scl = (-1)^(k*(k-1)/2)*a.scl
	return b
end
#########################################################
function bdual(a::kbase)
	l = length(a.e)
	e = Vector{Bool}(l)
	fill!(e, true)
	I = kbase(e)
	Inv = breverse(I)
	dual = bgeoprod(a, Inv)
	return dual
end
#########################################################
function remove(A::kMultvec, b::Number)
	l = length(A.comp)
	aux = 0
	for i=1:l
		if A.comp[i].scl == b
			aux = aux+1
		end
	end
	B = Vector{kbase}(l-aux)
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
	AB = Vector{kbase}(l+m)
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
				AB[j] = bsum(aux, AB[i])
				aux = AB[j]
				AB[i].scl = 0.0
			end
		end
	end
	result = remove(kMultvec(AB), 0.0)
	if length(result.comp) != 0
		return result
	else
		return kMultvec(kbase(id, 0.0))
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
		B = kMultvec([kbase(id, false, 0.0)])
	end
	return mvsum(A,B)
end
#########################################################
function reduct(A::kMultvec)
	l = length(A.comp)
	if l != 0
		a = kMultvec([A.comp[1]])
		B = Vector{kbase}(l-1)
		for i=2:l
			B[i-1] = A.comp[i]
		end
		C = mvsum(a, kMultvec(B))
		return C
	else
		return kMultvec([kbase(id, 0.0)])
	end
end
#########################################################
function geoprod(A::kMultvec, B::kMultvec)
	l = length(A.comp)
	m = length(B.comp)
	AB = Matrix{kbase}(l,m)
	for i=1:l
		for j=1:m
			AB[i,j] = bgeoprod(A.comp[i], B.comp[j])
		end
	end
	p = kMultvec(reshape(AB,l*m))
	p = reduct(p)
	result = remove(p, 0.0)
	if length(result.comp) != 0
		return result
	else
		return kMultvec([kbase(id, 0.0)])
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
	AB = Matrix{kbase}(l,m)
	for i=1:l
		for j=1:m
			AB[i,j] = innprod(A.comp[i], B.comp[j])
		end
	end
	p = kMultvec(reshape(AB,l*m))
	p = reduct(p)
	result = remove(p, 0)
	if length(result.comp) != 0
		return result
	else
		return kMultvec([kbase(id, 0.0)])
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
	AB = Matrix{kbase}(l,m)
	for i=1:l
		for j=1:m
			AB[i,j] = outprod(A.comp[i], B.comp[j])
		end
	end
	p = kMultvec(reshape(AB,l*m))
	p = reduct(p)
	result = remove(p, 0)
	if length(result.comp) != 0
		return result
	else
		return kMultvec([kbase(id, 0.0)])
	end
end
#########################################################
function Base.:^(A::kMultvec,B::kMultvec)
    return outer(A,B)
end
#########################################################
function copy(A::kMultvec)
	l = length(A.comp)
	X = Vector{kbase}(l)
	for i=1:l
		X[i] = copy(A.comp[i])
	end
	return kMultvec(X)
end
#########################################################
function dual(A::kMultvec)
	k = length(A.comp)
	dual = kMultvec(k)
	for i=1:k
		dual.comp[i] = bdual(A.comp[i])
	end
	return dual
end
#########################################################
function escalar(A::kMultvec, B::kMultvec)
	l = length(A.comp)
	m = length(B.comp)
	AB = Matrix{kbase}(l,m)
	for i=1:l
		for j=1:m
			AB[i,j] = escalarb(A.comp[i], B.comp[j])
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
    return escalar(A,B)
end
#########################################################
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
function escalar(A::kblade, B::kblade)
	X = bltomv(A)
	Y = bltomv(B)
	result = escalar(X, Y)
	return result
end
#########################################################
function Base.:*(A::kblade,B::kblade)
    return escalar(A,B)
end
#########################################################
function breverse(A::kblade)
	l = length(A.conj)
	T = Vector{kMultvec}(l)
	fill!(T, kMultvec(kbase(id, 0.0)))
	for i=1:l
		T[i] = A.conj[l-i+1]
	end
	return kblade(T)
end
#########################################################
function bnorm(A::kblade)
	X = escalar(A, breverse(A))
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
function inverse(A::kblade)
	div = escalar(A, breverse(A))
	l = length(A.conj)
	aux = breverse(A)
	X = copy(aux)
	m = length(aux.conj[1].comp)
	for j=1:m
		X.conj[1].comp[j].scl = (aux.conj[1].comp[j].scl)*(div)^(-1)
	end
	return X
end
#########################################################
function projection(A::kblade, N::kblade)
	N2 = inverse(N)
	result = geoprod(inner(A, N2), bltomv(N))
	return result
end
#########################################################
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
function show(io::IO, a::kbase)
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
type pbase
	eb::Vector{Bool}
	ep::Bool
	scl::Number

	pbase(a::Vector{Bool}, b::Bool, c::Number) = new(a, b, c)
	pbase(a::Vector{Bool}, b::Bool) = new(a, b, 1.0)
	pbase(b::Number) = new(id, false, b)
end
#########################################################
type pMultvec
	comp::Vector{pbase}

	pMultvec(a::pbase) = new([a])
	pMultvec(X::Vector{pbase}) = new(X)
	pMultvec(a::Number) = new([pbase(id, false, a)])
end
#########################################################
function grade(a::pbase)
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
function show(io::IO, a::pbase)
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
function bbtopb(a::kbase)
	b = pbase(a.e, false, a.scl)
	return b
end
#########################################################
function copy(a::pbase)
	eb = copy(a.eb)
	ep = copy(a.ep)
	esc = copy(a.scl)
	result = pbase(eb,ep, esc)
	return result
end
#########################################################
function bgeoprod(a::pbase, b::pbase)
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
	f = bgeoprod(kbase(c, a.scl), kbase(d, b.scl))
	scl = aux*f.scl
	g = Vector{Bool}(l)
	for i=1:l
		g[i] = f.e[i]
	end
	result = pbase(g, ep, scl)
	return result
end
#########################################################
function Base.:∘(a::pbase,b::pbase)
    return bgeoprod(A,B)
end
#########################################################
function innprod(a::pbase, b::pbase)
	if grade(bgeoprod(a,b)) == abs(grade(a) - grade(b))
		return bgeoprod(a,b)
	else
		x = Vector{Bool}(length(a.eb))
		fill!(x, false)
		return pbase(x, false, 0.0)
	end
end
#########################################################
function Base.:⋅(a::pbase,b::pbase)
    return innprod(A,B)
end
#########################################################
function outprod(a::pbase, b::pbase)
	if grade(bgeoprod(a,b)) == (grade(a) + grade(b))
		return bgeoprod(a,b)
	else
		x = Vector{Bool}(length(a.eb))
		fill!(x, false)
		return pbase(x, false, 0.0)
	end
end
#########################################################
function Base.:^(a::pbase,b::pbase)
    return outprod(A,B)
end
#########################################################
function escalarb(a::pbase, b::pbase)
	if grade(bgeoprod(a,b)) == 0
		return bgeoprod(a,b)
	else
		x = Vector{Bool}(length(a.eb))
		fill!(x, false)
		return pbase(x, false, 0.0)
	end
end
#########################################################
function Base.:*(a::pbase,b::pbase)
    return escalarb(A,B)
end
#########################################################
function bsum(a::pbase, b::pbase)
	l = length(a.eb)
	eb = Vector{Bool}(l)
	fill!(eb, false)
	if a.eb == b.eb && a.ep == b.ep
		scl = a.scl + b.scl
	 	return pbase(a.eb, a.ep, scl)
	else
		return pMultvec([a,b])
	end
end
#########################################################
function Base.:+(a::pbase,b::pbase)
    return bsum(a,b)
end
function Base.:-(a::pbase,b::pbase)
	b.scl = -b.scl
	return bsum(a,b)
end
#########################################################
function breverse(a::pbase)
	k = grade(a)
	eb = copy(a.eb)
	ep = copy(a.ep)
	b = pbase(eb,ep, 1)
	b.scl = (-1)^(k*(k-1)/2)*a.scl
	return b
end
#########################################################
function bconjugate(a::pbase)
	grm = 0
	if a.ep == true
		grm = 1
	end
	b = breverse(a)
	b.scl = b.scl*(-1)^grm
	return b
end
#########################################################
function bdual(a::pbase)
	l = length(a.eb)
	eb = Vector{Bool}(l)
	fill!(eb, true)
	I = pbase(eb,true,1.0)
	Inv = breverse(I)
	dual = bgeoprod(a, Inv)
	return dual
end
##################FUNÇÕES AUXILIARES#####################
function prtore(a::pbase)
	l = length(a.eb)
	e = Vector{Bool}(l+1)
	for i=1:l
		e[i] = a.eb[i]
	end
	e[l+1] = a.ep
	b = kbase(e, a.scl)
	return b
end
#########################################################
function retopr(a::kbase)
	l = length(a.e)
	eb = Vector{Bool}(l-1)
	ep = a.e[l]
	for i=1:l-1
		eb[i] = a.e[i]
	end
	b = pbase(eb,ep,a.scl)
end
#########################################################
function leneb(A::pMultvec)
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
function mvtopmv(A::kMultvec)
	l = length(A.comp)
	X = Vector{pbase}(l)
	for i=1:l
		X[i] = bbtopb(A.comp[i])
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
	B = Vector{pbase}(l-aux)
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
	X = Vector{pbase}(l)
	for i=1:l
		X[i] = copy(A.comp[i])
	end
	return pMultvec(X)
end
#########################################################
function mvsum(A::pMultvec, B::pMultvec)
	l = length(A.comp)
	m = length(B.comp)
	A2 = Vector{kbase}(l)
	B2 = Vector{kbase}(m)
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
	AB = Vector{pbase}(n)
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
		eb = leneb(B)
		B = pMultvec([pbase(eb, false, 0.0)])
	end
	return mvsum(A,B)
end
#######################AUXILIAR##########################
function reduct(A::pMultvec)
	l = length(A.comp)
	eb = leneb(A)
	if l != 0
		a = pMultvec([A.comp[1]])
		B = Vector{pbase}(l-1)
		for i=2:l
			B[i-1] = A.comp[i]
		end
		C = mvsum(a, pMultvec(B))
		return C
	else
		return pMultvec([pbase(eb, false, 0.0)])
	end
end
#########################################################
function geoprod(A::pMultvec, B::pMultvec)
	l = length(A.comp)
	m = length(B.comp)
	eb = leneb(A)
	AB = Matrix{pbase}(l,m)
	for i=1:l
		for j=1:m
			AB[i,j] = bgeoprod(A.comp[i], B.comp[j])
		end
	end
	p = pMultvec(reshape(AB,l*m))
	p = reduct(p)
	result = remove(p, 0)
	if length(result.comp) != 0
		return result
	else
		return pMultvec([pbase(eb, false, 0.0)])
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
	eb = leneb(A)
	AB = Matrix{pbase}(l,m)
	for i=1:l
		for j=1:m
			AB[i,j] = innprod(A.comp[i], B.comp[j])
		end
	end
	p = pMultvec(reshape(AB,l*m))
	p = reduct(p)
	result = remove(p, 0)
	if length(result.comp) != 0
		return result
	else
		return pMultvec([pbase(eb, false, 0.0)])
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
	AB = Matrix{pbase}(l,m)
	eb = leneb(A)
	for i=1:l
		for j=1:m
			AB[i,j] = outprod(A.comp[i], B.comp[j])
		end
	end
	p = pMultvec(reshape(AB,l*m))
	p = reduct(p)
	result = remove(p, 0)
	if length(result.comp) != 0
		return result
	else
		return pMultvec([pbase(eb, false, 0.0)])
	end
end
#########################################################
function Base.:^(A::pMultvec,B::pMultvec)
    return outer(A,B)
end
#########################################################
function escalar(A::pMultvec, B::pMultvec)
	l = length(A.comp)
	m = length(B.comp)
	AB = Matrix{pbase}(l,m)
	for i=1:l
		for j=1:m
			AB[i,j] = escalarb(A.comp[i], B.comp[j])
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
    return escalar(A,B)
end
#########################################################
function conjugate(A::pMultvec)
	l = length(A.comp)
	B = Vector{pbase}(l)
	for i=1:l
		B[i] = bconjugate(A.comp[i])
	end
	return pMultvec(B)
end
#########################################################
function mvnorm(A::pMultvec)
	scl = escalar(A,conjugate(A))
	result = sqrt(scl)
	return result
end
#########################################################
function dual(A::pMultvec)
	k = length(A.comp)
	dual = Vector{pbase}(k)
	for i=1:k
		dual[i] = bdual(A.comp[i])
	end
	return pMultvec(dual)
end
#########################################################
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
function escalar(A::pblade, B::pblade)
	X = bltomv(A)
	Y = bltomv(B)
	result = escalar(X, Y)
	return result
end
#########################################################
function Base.:*(A::pblade, B::pblade)
    return escalar(A,B)
end
#########################################################
function leneb(A::pblade)
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
function breverse(A::pblade)
	l = length(A.conj)
	eb = leneb(A)
	T = Vector{pMultvec}(l)
	fill!(T, pMultvec([pbase(eb, false, 0.0)]))
	for i=1:l
		T[i] = A.conj[l-i+1]
	end
	return pblade(T)
end
#########################################################
function conjugate(A::pblade)
	l = length(A.conj)
	T = Vector{pMultvec}(l)
	eb = leneb(A)
	fill!(T, pMultvec([pbase(eb, false, 0.0)]))
	for i=1:l
		T[i] = conjugate(A.conj[l-i+1])
	end
	return pblade(T)
end
#########################################################
function bnorm(A::pblade)
	X = escalar(A, conjugate(A))
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
		div = escalar(A, breverse(A))
		l = length(A.conj)
		aux = breverse(A)
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
type cbase
	er::Vector{Bool}
	ei::Bool
	eo::Bool
	scl::Number
	function cbase(a::Vector{Bool},b::Bool,c::Bool,d::Number)
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
type cMultvec
	comp::Vector{cbase}

end
#########################################################
function grade(a::cbase)
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
function show(io::IO, a::cbase)
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
	X = Vector{cbase}(l)
	for i=1:l
		X[i] = copy(A.comp[i])
	end
	return cMultvec(X)
end
#########################################################
function bgeoprod(a::cbase, b::cbase)
	l=length(a.er)
	aux = 1
	a1 = kbase(a.er, a.scl)
	b1 = kbase(b.er, b.scl)
	ab1 = bgeoprod(a1, b1)
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
	ab = cbase(ab1.e, ei, eo, aux*ab1.scl)
	return ab
end
#########################################################
function Base.:∘(a::cbase, b::cbase)
    return bgeoprod(a,b)
end
#########################################################
function inner(a::cbase, b::cbase)
	if grade(bgeoprod(a,b)) == abs(grade(a) - grade(b))
		return bgeoprod(a,b)
	else
		return cbase(id, false, false, 0.0)
	end
end
#########################################################
function Base.:⋅(a::cbase, b::cbase)
    return inner(a,b)
end
#########################################################
function outer(a::cbase, b::cbase)
	if grade(bgeoprod(a,b)) == grade(a) + grade(b)
		return bgeoprod(a,b)
	else
		return cbase(id, false, false, 0.0)
	end
end
#########################################################
function Base.:^(a::cbase, b::cbase)
    return outer(a,b)
end
#########################################################
function escalar(a::cbase, b::cbase)
	if grade(bgeoprod(a,b)) == 0
		return bgeoprod(a,b)
	else
		return cbase(id, false, false, 0.0)
	end
end
#########################################################
function Base.:*(a::cbase, b::cbase)
    return escalar(a,b)
end
#########################################################
function bsum(a::cbase, b::cbase)
	l = length(a.er)
	er = Vector{Bool}(l)
	fill!(er, false)
	if a.er == b.er && a.ei == b.ei && a.eo == b.eo
		scl = a.scl + b.scl
	 	return cbase(a.er, a.ei,a.eo, scl)
	else
		return cMultvec([a,b])
	end
end
#########################################################
function Base.:+(a::cbase, b::cbase)
    return bsum(a,b)
end
function Base.:-(a::cbase, b::cbase)
	b.scl = - b.scl
    return bsum(a,b)
end
#########################################################
function cbtopb(a::cbase)
	l = length(a.er)
	ea = fill!(Vector{Bool}(l+1), false)
	eb = copy(ea)
	eb[l+1] = true
	X = pMultvec(pbase(ea, false, 0.0))
	Y = copy(X)
	f = 0
	er = Vector{Bool}(l+1)
	for i=1:l
		er[i] = a.er[i]
	end
	er[l+1] = false
	if a.ei == true
		X = pMultvec([pbase(ea, true, 1.0), pbase(eb, false, 1.0)])
		f = 1
	end
	if a.eo == true
		Y = pMultvec([pbase(ea, true, 0.5), pbase(eb, false, -0.5)])
		f = 1
	end
	if f == 1
		p1 = geoprod(pMultvec(pbase(er, false, a.scl)), X)
		p2 = geoprod(pMultvec(pbase(er, false, a.scl)), Y)
		return mvsum(p1, p2)
	else
		return pbase(er, false, a.scl)
	end
end
#########################################################
function cbtore(a::cbase)
	l = length(a.er)
	e = Vector{Bool}(l+2)
	for i=1:l
		e[i] = a.er[i]
	end
	e[l+1] = a.ei
	e[l+2] = a.eo
	b = kbase(e, a.scl)
	return b
end
#########################################################
function retocb(a::kbase)
	l = length(a.e)
	er = Vector{Bool}(l-2)
	for i=1:l-2
		er[i] = a.e[i]
	end
	ei = a.e[l-1]
	eo = a.e[l]
	b = cbase(er, ei, eo, a.scl)
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
	B = Vector{cbase}(l-aux)
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
	A2 = Vector{kbase}(l)
	B2 = Vector{kbase}(m)
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
	AB = Vector{cbase}(n)
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
		B = cMultvec([cbase(id,false, false, 0.0)])
	end
	return mvsum(A,B)
end
#########################################################
function reduct(A::cMultvec)
	l = length(A.comp)
	if l != 0
		a = cMultvec([A.comp[1]])
		B = Vector{cbase}(l-1)
		for i=2:l
			B[i-1] = A.comp[i]
		end
		C = mvsum(a, cMultvec(B))
		return C
	else
		return cMultvec([cbase(id,false, false, 0.0)])
	end
end
#########################################################
function geoprod(A::cMultvec, B::cMultvec)
	l = length(A.comp)
	m = length(B.comp)
	AB = Matrix{cbase}(l,m)
	for i=1:l
		for j=1:m
			AB[i,j] = bgeoprod(A.comp[i], B.comp[j])
		end
	end
	p = cMultvec(reshape(AB,l*m))
	p = reduct(p)
	result = remove(p, 0)
	if length(result.comp) != 0
		return result
	else
		return cMultvec([cbase(id,false, false, 0.0)])
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
	AB = Matrix{cbase}(l,m)
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
		return cMultvec([cbase(id,false, false, 0.0)])
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
	AB = Matrix{cbase}(l,m)
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
		return cMultvec([cbase(id,false, false, 0.0)])
	end
end
#########################################################
function Base.:^(A::cMultvec,B::cMultvec)
    return outer(A,B)
end
#########################################################
function escalar(A::cMultvec, B::cMultvec)
	l = length(A.comp)
	m = length(B.comp)
	AB = Matrix{cbase}(l,m)
	for i=1:l
		for j=1:m
			AB[i,j] = escalar(A.comp[i], B.comp[j])
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
    return escalar(A,B)
end
#########################################################
function pbtocb(a::pbase)
	l = length(a.eb)
	eb = fill!(Vector{Bool}(l-1), false)
	for i=1:l-1
		eb[i] = a.eb[i]
	end
	a2 = cbase(eb, false, false, a.scl)
	aux1 = cMultvec([cbase(id, false, false, 1.0)])
	aux2 = cMultvec([cbase(id, false, false, 1.0)])
	f = 0
	if a.eb[l] == true
		aux1 = cMultvec([cbase(id, true, false, 0.5), cbase(id,false, true, -1.0)])
		f = 1
	end
	if a.ep == true
		aux2 = cMultvec([cbase(id, true, false, 0.5), cbase(id, false, true, 1.0)])
		f = 1
	end
	b = cbase(eb, false, false, a.scl)
	if f == 1
		b = cMultvec([b])
		p1 = geoprod(aux1, aux2)
		return geoprod(b, p1)
	else
		return b
	end
end
#########################################################
function copy(a::cbase)
	er = copy(a.er)
	ei = copy(a.ei)
	eo = copy(a.eo)
	scl = copy(a.scl)
	return cbase(er,ei,eo,scl)
end
#########################################################
function pvtocv(A::pMultvec)
	l = length(A.comp)
	B = Vector{Any}(l)
	for i=1:l
		aux = pbtocb(A.comp[i])
		if typeof(aux) == cbase
			B[i] = aux
		else
			B[i] = cbase(id, false,false, 0.0)
			B = vcat(B, aux.comp)
		end
	end
	B = cMultvec(B)
	B = reduct(B)
	return B
end
#########################################################
function cvtopv(A::cMultvec)
	l = length(A.comp)
	m = length(A.comp[1].er)
	eb = fill!(Vector{Bool}(m+1), false)
	B = Vector{pbase}(l)
	for i=1:l
		aux = cbtopb(A.comp[i])
		if typeof(aux) == pbase
			B[i] = aux
		else
			B[i] = pbase(eb,false, 0.0)
			B = vcat(B, aux.comp)
		end
	end
	B = pMultvec(B)
	B = reduct(B)
	return pMultvec(B)
end
#########################################################
function retoaffin(a::Vector{Float64})
	l = length(a)
	A = Vector{kbase}(l+1)
	aux = Vector{Bool}(l+1)
	fill!(aux, false)
	for i=1:l
		aux[i] = true
		A[i] = copy(kbase(aux, a[i]))
		aux[i] = false
	end
	aux[l+1] = true
	A[l+1] = copy(kbase(aux, 1.0))
	return kMultvec(A)
end
#########################################################
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
		scl = inner(A, kMultvec(kbase(aux, 1.0))).comp[1].scl
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
		X.comp[aux] = kbase(e, 0.0)
		return reduct(X)
	elseif f1 == 1 && f2 == 0
		error("Image error")
	elseif f1 == 0 && f2 == 1
		error("Grade error")
	end
end
#########################################################
function euctoga(x::Vector{Float64})
	l = length(x)
	X = Vector{kbase}(l)
	e = fill!(Vector{Bool}(l), false)
	for i=1:l
		e[i] = true
		X[i] = kbase(copy(e), x[i])
		e[i] = false
	end
	return kMultvec(X)
end
#########################################################
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
function H(x::Vector{Float64})
	l = length(x)
	y = euctoga(x)
	X = Vector{pbase}(l+1)
	for i=1:l
		X[i] = bbtopb(y.comp[i])
	end
	m = length(S.comp[1].e)
	ide = Vector{Bool}(m+1)
	fill!(ide, false)
	X[l+1] = pbase(ide, true, 1.0)
	return pMultvec(X)
end
#########################################################
function H(S::kMultvec)
	l = length(S.comp)
	X = Vector{pbase}(l+1)
	for i=1:l
		X[i] = bbtopb(S.comp[i])
	end
	m = length(S.comp[1].e)
	ide = Vector{Bool}(m)
	fill!(ide, false)
	X[l+1] = pbase(ide, true, 1.0)
	return pMultvec(X)
end
#########################################################
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
#########################################################
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
	ei = pMultvec([pbase(eb, true, 1.0), pbase(eb2, false, 1.0)])
	eo = pMultvec([pbase(eb, true, 0.5), pbase(eb2, false, -0.5)])
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
function iconformal(X::cMultvec)
	Y = cvtopv(X)
	Y = ipconformal(Y)
	return pvtocv(Y)
end
#########################################################
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
function bltomv(A::cblade)
	l = length(A.conj)
	X = A.conj[1]
	for i=2:l
		X = outer(X, A.conj[i])
	end
	return X
end
#########################################################
function Base.:∘(A::cblade,B::cblade)
	AA = bltomv(A)
	BB = bltomv(B)
	return geoprod(AA,BB)
end
#########################################################
function Base.:⋅(A::cblade,B::cblade)
	AA = bltomv(A)
	BB = bltomv(B)
	return inner(AA,BB)
end
#########################################################
function Base.:^(A::cblade,B::cblade)
	AA = bltomv(A)
	BB = bltomv(B)
	return outer(AA,BB)
end
#########################################################
function Base.:*(A::cblade,B::cblade)
	AA = bltomv(A)
	BB = bltomv(B)
	return escalar(AA,BB)
end
#########################################################
function Base.:+(A::cblade,B::cblade)
	AA = bltomv(A)
	BB = bltomv(B)
	return mvsum(AA,BB)
end
#########################################################
function Base.:-(A::cblade,B::cblade)
	AA = bltomv(A)
	BB = bltomv(B)
	return AA-BB
end
#########################################################
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
function breverse(A::cblade)
	B = cbltopbl(A)
	B = breverse(B)
	return pbltocbl(B)
end
#########################################################
function conjugate(A::cblade)
	B = cbltopbl(A)
	B = conjugate(B)
	return pbltocbl(B)
end
#########################################################
function bnorm(A::cblade)
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
		div = (A * breverse(A))
		l = length(A.conj)
		aux = breverse(A)
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
