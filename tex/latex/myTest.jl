using Base
using DataFrames
using Distributions
using Gadfly
import Base




type MyObject
  x::Integer
  y::Integer
end
a=MyObject(23,44)

+(a::MyObject,b::MyObject)= MyObject(a.x+b.x,a.y,b.x)

generic size(obj::Any, [dim::Integer])

size(a::Array) = arraysize(a)


size(a::Array, d) = arraysize(a, d)
