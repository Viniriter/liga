import Base.show, Base.copy

type bblade
    e::Vector{Bool}
	escalar::Number

	bblade(a::Vector{Bool}) =   new(a, 1.0)
	bblade(a::Vector{Bool}, b::Number) =   new(a, b)
end

type MultVec
    comp::Vector{bblade}

	MultVec(a::bblade) = new([a])
	MultVec(A::Vector{bblade}) = new(A)
	MultVec(n::Int) = new(Vector{bblade}(n))
end

#######GRADE DOS ELEMENTOS BÁSICOS#######
function grade(x::bblade)
	acu = 0
	for i=1:length(x.e)
		if x.e[i] == true
			acu = acu + 1
		end
	end
	return acu
end
#########################################
function mvectovec(X::Vector{MultVec})
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
						A[i,k] = X[i].comp[j].escalar
					end
				end
			end
		end
		return A
	end
end

type blade
	conj::Vector{MultVec}

	function blade(X::Vector{MultVec})
		T = mvectovec(X)
		if det(T*transpose(T)) != 0
			new(X)
		else
			error("L.d")
		end
	end
end

#############VIZUALIZAÇÃO################
function show(io::IO, a::bblade)
	if (a.escalar == -1) && (grade(a) == 0)
		print(io, "-1")
	elseif (a.escalar == -1) && (grade(a) != 0)
		print(io, "-")
	elseif (a.escalar != 1) || (grade(a) == 0)
		print(io, "$(a.escalar)")
	end
	flag = 0
	for i=1:length(a.e)
		if a.e[i] == true
			flag = 1
		end
	end
	if flag == 1 && a.escalar != 0
		print(io, "e")
	end
	for i=1:length(a.e)
		if a.e[i] == true && a.escalar != 0
			print(io, "$i")
		end
	end
end

function show(io::IO, A::MultVec)
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

function show(io::IO, X::blade)
	l = length(X.conj)
	for i=1:l-1
		print(io, "($(X.conj[i]))∧")
	end
	print(io, "($(X.conj[l]))")
end
#########################################
######OPERAÇÕES COM BLADES BÁSICAS#######
#########################################
function copy(a::bblade)
	e = copy(a.e)
	esc = copy(a.escalar)
	result = bblade(e, esc)
	return result
end

#######PRODUTO GEOMÉTRICO BASE###########
function bgeoprod(a::bblade, b::bblade)
	l = length(a.e)
	Ea = copy(a.e)
	Eb = copy(b.e)
	ab = bblade(copy(Ea), 0)
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
	ab.escalar = a.escalar*b.escalar*(-1)^cont
	return ab
end
#########################################

#########PRODUTO INTERNO BASE############
function innprod(a::bblade, b::bblade)
	if grade(bgeoprod(a,b)) == abs(grade(a) - grade(b))
		return bgeoprod(a,b)
	else
		x = Vector{Bool}(length(a.e))
		fill!(x, false)
		return bblade(x, 0)
	end
end
#########################################

########PRODUTO INTERNO BASE#############
function outprod(a::bblade, b::bblade)
	if grade(bgeoprod(a,b)) == (grade(a) + grade(b))
		return bgeoprod(a,b)
	else
		x = Vector{Bool}(length(a.e))
		fill!(x, false)
		return bblade(x, 0)
	end
end
#########################################
function escalarb(a::bblade, b::bblade)
	if grade(bgeoprod(a,b)) == 0
		return bgeoprod(a,b)
	else
		x = Vector{Bool}(length(a.e))
		fill!(x, false)
		return bblade(x, 0)
	end
end
#################SOMA####################
function bsum(a::bblade, b::bblade)
	l = length(a.e)
	e = Vector{Bool}(l)
	fill!(e, false)
	sume = bblade(e, 0)
	if a.e == b.e
		sume.e = a.e
		sume.escalar = a.escalar + b.escalar
	elseif a.e != b.e
		sume = MultVec([a,b])
	end
	return sume
end
#########################################

################REVERSO##################
function breverse(a::bblade)
	k = grade(a)
	e = copy(a.e)
	b = bblade(e)
	b.escalar = (-1)^(k*(k-1)/2)*a.escalar
	return b
end
#########################################

###############CONJUGADO#################
function conjugate(a::bblade)
 #só faz sendido em G(p,q)
end

#########################################
function bdual(a::bblade)
	l = length(a.e)
	e = Vector{Bool}(l)
	fill!(e, true)
	I = bblade(e)
	Inv = breverse(I)
	dual = bgeoprod(a, Inv)
	return dual
end
#########################################

#########################################
######$OPERAÇÕES COM MULTIVETORES$#######
#########################################

############SOMA MULTIVETORES############
function mvsum(A::MultVec, B::MultVec)
	l = length(A.comp)
	m = length(B.comp)
	AB = Vector{bblade}(l+m)
	for i=1:l
		AB[i] = copy(A.comp[i])
	end
	for i=l+1:l+m
		AB[i] = copy(B.comp[i-l])
	end
	for j=1:(m+l-1)
		aux = AB[j]
		for i=j+1:m+l
			if AB[i].e == aux.e && AB[i].escalar != 0.0
				AB[j] = bsum(aux, AB[i])
				AB[i].escalar = 0.0
			end
		end
	end
	return MultVec(AB)
end
#########################################

###########PRODUTO GEOMÉTRICO############
function geoprod(A::MultVec, B::MultVec)
	l = length(A.comp)
	m = length(B.comp)
	acu = Matrix(l,m)
	for i=1:l
		for j=1:m
			acu[i,j] = bgeoprod(A.comp[i],B.comp[j])
		end
	end
	p = reshape(acu, l*m)
	p2 = Vector{bblade}(l*m)
	for i=1:l*m
		p2[i] = p[i]
	end
	p = MultVec(p2)
	rr = MultVec(p.comp[1])
	for i=2:(l*m)
		rr = mvsum(rr, MultVec(p.comp[i]))
	end
	zeros = 0
	for i=1:l*m
		if rr.comp[i].escalar == 0
			zeros = zeros + 1
		end
	end
	result = MultVec(l*m - zeros)
	for i=1:(l*m)
		for j=1:l*m
			if abs(rr.comp[i].escalar) >= abs(rr.comp[j].escalar)
				aux = rr.comp[j]
				rr.comp[j] = rr.comp[i]
				rr.comp[i] = aux
			end
		end
	end
	for i=1:(l*m-zeros)
		result.comp[i] = rr.comp[i]
	end
	for i=1:(l*m-zeros)
		for j=1:(l*m-zeros)
			if grade(result.comp[i]) <= grade(result.comp[j])
				aux = result.comp[i]
				result.comp[i] = result.comp[j]
				result.comp[j] = aux
			end
		end
	end
	return result
end
#########################################

############PRODUTO INTERNO##############
function inner(A::MultVec, B::MultVec)
	l = length(A.comp)
	m = length(B.comp)
	acu = Matrix(l,m)
	for i=1:l
		for j=1:m
			acu[i,j] = innprod(A.comp[i],B.comp[j])
		end
	end
	p = reshape(acu, l*m)
	p2 = Vector{bblade}(l*m)
	for i=1:l*m
		p2[i] = p[i]
	end
	p = MultVec(p2)
	p2 = MultVec(copy(p.comp))
	rr = MultVec(p2.comp[1])
	for i=2:(l*m)
		rr = mvsum(rr, MultVec(p2.comp[i]))
	end
	for i=1:(l*m)
		for j=1:l*m
			if abs(rr.comp[i].escalar) >= abs(rr.comp[j].escalar)
				aux = rr.comp[j]
				rr.comp[j] = rr.comp[i]
				rr.comp[i] = aux
			end
		end
	end
	zeros = 0
	for i=1:l*m
		if rr.comp[i].escalar == 0
			zeros = zeros + 1
		end
	end
	result = MultVec(l*m-zeros)
	for i=1:(l*m-zeros)
		result.comp[i] = rr.comp[i]
	end
	return result
end
#########################################

############PRODUTO EXTERIOR#############
function outer(A::MultVec, B::MultVec)
	l = length(A.comp)
	m = length(B.comp)
	acu = Matrix(l,m)
	for i=1:l
		for j=1:m
			acu[i,j] = outprod(A.comp[i],B.comp[j])
		end
	end
	p = reshape(acu, l*m)
	p2 = Vector{bblade}(l*m)
	for i=1:l*m
		p2[i] = p[i]
	end
	p = MultVec(p2)
	p2 = MultVec(copy(p.comp))
	rr = MultVec(p2.comp[1])
	for i=2:(l*m)
		rr = mvsum(rr, MultVec(p2.comp[i]))
	end
	for i=1:(l*m)
		for j=1:l*m
			if abs(rr.comp[i].escalar) >= abs(rr.comp[j].escalar)
				aux = rr.comp[j]
				rr.comp[j] = rr.comp[i]
				rr.comp[i] = aux
			end
		end
	end
	zeros = 0
	for i=1:l*m
		if rr.comp[i].escalar == 0
			zeros = zeros + 1
		end
	end
	result = MultVec(l*m-zeros)
	for i=1:(l*m-zeros)
		result.comp[i] = rr.comp[i]
	end
	return result
end
#########################################
function copy(A::MultVec)
	l = length(A.comp)
	X = Vector{bblade}(l)
	for i=1:l
		X[i] = copy(A.comp[i])
	end
	return MultVec(X)
end
##################DUAL###################
function dual(A::MultVec)
	k = length(A.comp)
	dual = MultVec(k)
	for i=1:k
		dual.comp[i] = bdual(A.comp[i])
	end
	return dual
end
#########################################
function escalar(A::MultVec, B::MultVec)
	l = length(A.comp)
	m = length(B.comp)
	acu = Matrix(l,m)
	for i=1:l
		for j=1:m
			acu[i,j] = escalarb(A.comp[i],B.comp[j])
		end
	end
	p = reshape(acu, l*m)
	p2 = Vector{bblade}(l*m)
	for i=1:l*m
		p2[i] = p[i]
	end
	p = MultVec(p2)
	p2 = MultVec(copy(p.comp))
	rr = MultVec(p2.comp[1])
	for i=2:(l*m)
		rr = mvsum(rr, MultVec(p2.comp[i]))
	end
	for i=1:(l*m)
		for j=1:l*m
			if abs(rr.comp[i].escalar) >= abs(rr.comp[j].escalar)
				aux = rr.comp[j]
				rr.comp[j] = rr.comp[i]
				rr.comp[i] = aux
			end
		end
	end
	zeros = 0
	for i=1:l*m
		if rr.comp[i].escalar == 0
			zeros = zeros + 1
		end
	end
	result = MultVec(l*m-zeros)
	for i=1:(l*m-zeros)
		result.comp[i] = rr.comp[i]
	end
	if length(result.comp) != 0
		return result.comp[1].escalar
	else
		return 0
	end
end
################NORMA####################
function mvnorm(A::MultVec)
	scl = escalar(A,A)
	result = sqrt(scl^2)
	return result
end
#########################################
#########################################
################BLADES###################
#########################################
function bltomv(A::blade)
	l = length(A.conj)
	X = A.conj[1]
	for i=2:l
		X = outer(X, A.conj[i])
	end
	return X
end
#########################################
function geoprod(A::blade, B::blade)
	X = bltomv(A)
	Y = bltomv(B)
	result = geoprod(X, Y)
	return result
end
#########################################
function outer(A::blade, B::blade)
	X = bltomv(A)
	Y = bltomv(B)
	result = outer(X, Y)
	return result
end
#########################################
function inner(A::blade, B::blade)
	X = bltomv(A)
	Y = bltomv(B)
	result = inner(X, Y)
	return result
end
#########################################
function escalar(A::blade, B::blade)
	X = bltomv(A)
	Y = bltomv(B)
	result = escalar(X, Y)
	return result
end
#########################################
function breverse(A::blade)
	l = length(A.conj)
	T = Vector{MultVec}(l)
	fill!(T, MultVec(bblade(id)))
	for i=1:l
		T[i] = A.conj[l-i+1]
	end
	return blade(T)
end
#########################################
function bnorm(A::blade)
	X = escalar(A, breverse(A))
	return sqrt(X)
end
#########################################
function dual(A::blade)
	X = bltomv(A)
	return dual(X)
end
#########################################
function copy(A::blade)
	l = length(A.conj)
	X = Vector{MultVec}(l)
	for i=1:l
		X[i] = copy(A.conj[i])
	end
	return blade(X)
end
#########################################
################INVERSO##################
#########################################
function inverse(A::blade)
	div = geoprod(A, breverse(A)).comp[1].escalar
	l = length(A.conj)
	aux = breverse(A)
	X = copy(aux)
	m = length(aux.conj[1].comp)
	for j=1:m
		X.conj[1].comp[j].escalar = (aux.conj[1].comp[j].escalar)*(div)^(-1)
	end
	return X
end
#########################################



#########################################################################
#########################################################################
################################TESTES###################################
#########################################################################
#########################################################################
include("layout3.jl")
layout(3)


a = MultVec([bblade(e1,2.0),bblade(e3,4.0)])
b = MultVec([bblade(e1,1.0), bblade(e2), bblade(e3,2.0)])
c = MultVec([bblade(e1)])
d = MultVec([bblade(e2)])
f = MultVec([bblade(e3)])
A = blade([a,b])
B = blade([c,d,f])
bltomv(A)
bltomv(breverse(A))
bltomv(B)
geoprod(A,B)
outer(A,B)
inner(A,B)
escalar(A,B)
breverse(A)
breverse(B)
geoprod(A, breverse(A))
bnorm(A)
bnorm(B)
inverse(A)
geoprod(A, inverse(A))
geoprod(inverse(A), A)
