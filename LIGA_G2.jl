importall Base.Operators

const basis = Dict("id"=>1,"e1"=>2, "e2"=>3, "e12"=>4)

const elb=["id","e1","e2","e12"]

const n=2

abstract type PairType end
struct PairT <: PairType
    sig :: Number
    val :: String
end

abstract type MultiVectorG2Type end
struct MVecG2 <: MultiVectorG2Type
    id :: Number
    e1 :: Number
    e2 :: Number
    e12:: Number
end

function basisprod(elbasisf,elbasiss)
    bp=[PairT(1,"id")  PairT(1,"e1")   PairT(1,"e2")   PairT(1,"e12");
        PairT(1,"e1")  PairT(1,"id")   PairT(1,"e12")  PairT(1,"e2");
        PairT(1,"e2")  PairT(-1,"e12") PairT(1,"id")   PairT(-1,"e1");
        PairT(1,"e12") PairT(-1,"e2")  PairT(1,"e1")   PairT(-1,"id")]
    return bp[basis[elbasisf],basis[elbasiss]]
end 

function gaprod(u::MVecG2,v::MVecG2)
    uext=[u.id u.e1 u.e2 u.e12]
    vext=[v.id v.e1 v.e2 v.e12]   
    posi=0
    posj=0
    dimprob=2^n
    gp=Array{PairT}(dimprob*dimprob)
    posgp=1
    for i in uext
        posi+=1
        posj=0
        for j in vext
            posj+=1
            gp[posgp]=PairT(i*j*basisprod(elb[posi],elb[posj]).sig,basisprod(elb[posi],elb[posj]).val)
            posgp+=1
        end
    end    
    mvecres=zeros(dimprob)
    for i=1:dimprob*dimprob
        for j=1:dimprob
            if gp[i].val==elb[j]
             mvecres[j]=mvecres[j]+gp[i].sig
            end
        end
    end
    return MVecG2(mvecres[1],mvecres[2],mvecres[3],mvecres[4])
end

function Base.:+(u::MVecG2,v::MVecG2)
    return MVecG2(u.id+v.id,u.e1+v.e1,u.e2+v.e2,u.e12+v.e12)
end

function Base.:-(u::MVecG2,v::MVecG2)
    return MVecG2(u.id-v.id,u.e1-v.e1,u.e2-v.e2,u.e12-v.e12)
end

function Base.:∘(u::MVecG2,v::MVecG2)
    return gaprod(u::MVecG2,v::MVecG2)
end

function Base.:*(a::Number,u::MVecG2)
    return MVecG2(a*(u.id),a*(u.e1),a*(u.e2),a*(u.e12))
end

function Base.:*(u::MVecG2,a::Number)
    return MVecG2(a*u.id,a*u.e1,a*u.e2,a*u.e12)
end

function escprod(u::MVecG2,v::MVecG2)
    return (0.5)*(gaprod(u,v)+gaprod(v,u))
end

function outprod(u::MVecG2,v::MVecG2)
    return (0.5)*(gaprod(u,v)-gaprod(v,u))
end

function Base.:*(u::MVecG2,v::MVecG2)
    return escprod(u,v)
end

function Base.:^(u::MVecG2,v::MVecG2)
    return outprod(u,v)
end

import Base.show
function Base.show(io::IO, MV::MVecG2)
    print(io, "($(MV.id))id +($(MV.e1)) e1+ ($(MV.e2)) e2 +($(MV.e12)) e12")
end