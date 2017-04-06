using Gadfly, Distributions,DataArrays,DataFrames, ColorBrewer, Cairo

generation = 500 #  number of generations
rep_nb= 100      # number of repetitions
X=[]
matrix_mutation= zeros(100,500)
generation_vector = zeros(rep_nb,1) #timevector t=(x1,x2,...,xN)
matrix_pure= zeros(rep_nb,generation) #evolution matrix

#create dataframe for statistics
ev_data = DataFrame(Any[Pop_size=Float64[], mean=Float64[], std=Float64[]])

#plotting table
blankTheme = Theme(
grid_color=colorant"grey",
panel_fill=colorant"white",
major_label_color=colorant"black",
minor_label_color=colorant"black")

#simulate pure genetic drift
function gen_drift_p(a::Matrix, N::Integer, number::Integer)
  a[:,1] = N/2
  for rep in 1:rep_nb
    for gen in 2:generation
      #next generation will have random number of individuals A with probability of current generation
      p= matrix_pure[rep,gen-1]/N
      A= rand(Binomial(N,p))

      #update state
      a[rep,gen]=A

      #save the generation number and break if extinction/fixation occurs
      if a[rep,gen]==0 || matrix_pure[rep,gen]==N
        generation_vector[rep]= gen
        a[rep,gen+1:generation]= NaN
        break
      end
    end
  end
  #plotting window (module @layers using ColorBrewer using Gadfly end...???)
  colors = palette("Spectral", 11)
  lay =[layer(x=1:generation, y=matrix_pure[i,:], Geom.path,Theme(default_color=colors[i])) for i in 1:number] #Guide.manual_color_key("legend for plot", ["1","2"]["deepskyblue", "orange"])

  #plot
  plt=plot(lay... , Guide.title("genetic drift"),  Guide.XLabel("generation number"),Guide.YLabel("population number"),blankTheme )
  img = SVG("image/GenDrift_pure/population_size_$N.svg", 8inch, 6inch)
  draw(img, plt)
  return a
end
gen_drift_p(test,100,10)
 #evaluate data
function get_data(f::Function, df::DataFrame, N::Float64)
  if nrow(df) >0
    deleterows!(df, 1:nrow(df))
  end
    for i in 100:100:N
      f(i,10)
      push!(df, [i floor(mean(generation_vector)) std(generation_vector)])
      #push!(X, floor(mean(generation_vector)))
    end
    return df
end

function gen_drift_m(m::Matrix, N::Integer,number::Integer)
  μ=0.001 #mutationrate
  m[:,1] = N/2
  for rep in 1:rep_nb
    for gen in 2:generation

      #add mutation
      a=m[rep,gen-1]
      a_p=rand(Poisson(μ*a))
      b=N-m[rep,gen-1]
      b_p=rand(Poisson(μ*b))
      mutate = b_p-a_p
      m[rep,gen]=m[rep,gen-1]+mutate

      #genetic drift loop
      p= m[rep,gen]/N
      A= rand(Binomial(N,p))

      #update state
      m[rep,gen]=A

      #count extinction/fixation
      if m[rep,gen]==0 || m[rep,gen]==N
        generation_vector[rep]= gen

      end
    end
  end
  #plotting window (module @layers using ColorBrewer using Gadfly end...???)
  colors = palette("Spectral", 11)
  lay =[layer(x=1:generation, y=m[i,:], Geom.path,Theme(default_color=colors[i])) for i in 1:number] #Guide.manual_color_key("legend for plot", ["1","2"]["deepskyblue", "orange"])
  #plot
  plt=plot(lay... , Guide.title("genetic drift with mutationrate μ = $μ "),  Guide.XLabel("generation number"),Guide.YLabel("population size"),blankTheme )
  img = SVG("image/GenDrift_mutation/$N mutation_rate_$μ .svg", 8inch, 6inch)
  draw(img, plt)
  return m
end

get_data(gen_drift_m,ev_data,1000)
plot(ev_data, x="x2", Geom.density, Guide.title("Boxplot of mean value fixation occurs"), Guide.XLabel("population size"), blankTheme)
f=convert(Array{Integer}, generation_vector)
plot(matrix_pure, x=matrix_pure[:,500], Geom.histogram(bincount=30, density=true), blankTheme)
gen_drift_m(100,10)
