---
layout: page
title: Getting Started
navigation: 2
---

## Installation


(In the future) Open Julia REPL and type

    $ Pkg.add("liga")

Or (by now) install it yourself as:

    $ git clone https://github.com/evcastelani/liga

## Layouts

As cited, just the second option is avaiable. So, you can type 

	$ include("liga.jl")

This file contains all functions of the package. The first point that you must
be in mind is ```setup``` the layout. The ```layout.jl``` is a file that write
an ```objects.jl```. The ```object.jl``` generates a base for some space that 
we want work. For example, if we type

	$ layout(3)

we are generating  the G3 space, that is, we are generating a base of vectors and
multivectors for G3.
