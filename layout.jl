import Base.show

function layout(dim::Int)
    stf=readdlm("setup.dat")
    obj="objects.jl"
    f=open(obj,"w")
    println(f, "const n= $dim")
    println(f,"const bblade=$(stf[:,1])")
    println(f,"const bbval=$(stf[:,2])")
    mvtype="MultiVectorType"
    println(f,"abstract type $(mvtype) end")
    println(f,"struct MVec <: $(mvtype)")
    println(f,"coord::Array{Float64,1}")
    println(f,"end")
    println(f,"  ")
    println()
    println(f,"function Base.show(io::IO,mv::MVec)")
    println(f,"for i=1:2^n")
    println(f,"if mv.coord[i]>0")
    println(f,"   print(io,string( + ),mv.coord[i],bblade[i] )")
    println(f,"end")
    println(f,"if mv.coord[i]<0")
    println(f,"   print(io,string( ),mv.coord[i],bblade[i] )")
    println(f,"end")
    println(f,"end")
    println(f,"end")
    close(f)
    include(obj)
end