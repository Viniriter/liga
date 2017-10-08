---
layout: page
title: Getting Started
navigation: 2
---

## Installation


(In the future) Open Julia REPL and type

     Pkg.add("liga")

Or (by now) install it yourself as:

     git clone https://github.com/evcastelani/liga

## Layouts

As cited, just the second option is avaiable. So, you can type 

	 include("liga.jl")

This file contains all functions of the package. The first point that you must
be in mind is *setup* the layout. The ```layout.jl``` is a file that write
an ```objects.jl```. The ```objects.jl``` generates a basis for some space that we want work. For example, if we type

	 layout(3)

we are generating  the G3 space, that is, we are generating a basis of vectors and multivectors for G3. Now, was created bases elements in ```objects.jl```. Type, for example

     e1

This is a logical array with value 

     [false,true,true]

All bases elements were created using this idea. Try type others basis elements
    
     e123

     e12       
        
An  observation in this case is the order of these elements are important. If you type

     e21

an error mensage is showed. This does not mean that this element can not be built. We will see this later. Reciprocally, type

     kbase([true,true,false])

This command show the multivector ```e12```. Another important property is...    
