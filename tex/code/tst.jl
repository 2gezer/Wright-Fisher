using Gadfly, Distributions,DataArrays,DataFrames, ColorBrewer, Cairo

generation = 500 #  number of generations
rep_nb= 100      # number of repetitions
X=[]

generation_vector = zeros(rep_nb,1) #timevector t=(x1,x2,...,xN)
ev_matrix= zeros(rep_nb,generation) #evolution matrix

#create dataframe for statistics
ev_data = DataFrame(Any[Pop_size=Float64[], mean=Float64[], std=Float64[]])

#plotting table
blankTheme = Theme(
grid_color=colorant"grey",
panel_fill=colorant"white",
major_label_color=colorant"black",
minor_label_color=colorant"black")

#simulate pure genetic drift
function gen_drift_p(N,number)
  ev_matrix[:,1] = N/2
  for rep in 1:rep_nb
    for gen in 2:generation
      #next generation will have random number of individuals A with probability of current generation
      p= ev_matrix[rep,gen-1]/N
      A= rand(Binomial(N,p))

      #update state
      ev_matrix[rep,gen]=A

      #save the generation number and break if extinction/fixation occurs
      if ev_matrix[rep,gen]==0 || ev_matrix[rep,gen]==N
        generation_vector[rep]= gen
        ev_matrix[rep,gen+1:generation]= NaN
        break
      end
    end
  end
  #plotting window (module @layers using ColorBrewer using Gadfly end...???)
  colors = palette("Spectral", 11)
  lay =[layer(x=1:generation, y=ev_matrix[i,:], Geom.path,Theme(default_color=colors[i])) for i in 1:number] #Guide.manual_color_key("legend for plot", ["1","2"]["deepskyblue", "orange"])

  #plot
  plt=plot(lay... , Guide.title("genetic drift"),  Guide.XLabel("generation number"),Guide.YLabel("population number"),blankTheme )
  img = SVG("image/GenDrift_pure/population_size_$N.svg", 8inch, 6inch)
  draw(img, plt)
  return ev_matrix , ev_data
end

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

function gen_drift_m(N,number)
  μ=0.001 #mutationrate
  ev_matrix[:,1] = N/2
  for rep in 1:rep_nb
    for gen in 2:generation

      #add mutation
      a=ev_matrix[rep,gen-1]
      a_p=rand(Poisson(μ*a))
      b=N-ev_matrix[rep,gen-1]
      b_p=rand(Poisson(μ*b))
      mutate = b_p-a_p
      ev_matrix[rep,gen]=ev_matrix[rep,gen-1]+mutate

      #genetic drift loop
      p= ev_matrix[rep,gen]/N
      A= rand(Binomial(N,p))

      #update state
      ev_matrix[rep,gen]=A

      #count extinction/fixation
      if ev_matrix[rep,gen]==0 || ev_matrix[rep,gen]==N
        generation_vector[rep]= gen

      end
    end
  end
  #plotting window (module @layers using ColorBrewer using Gadfly end...???)
  colors = palette("Spectral", 11)
  lay =[layer(x=1:generation, y=ev_matrix[i,:], Geom.path,Theme(default_color=colors[i])) for i in 1:number] #Guide.manual_color_key("legend for plot", ["1","2"]["deepskyblue", "orange"])
  #plot
  plt=plot(lay... , Guide.title("genetic drift with mutationrate μ = $μ "),  Guide.XLabel("generation number"),Guide.YLabel("population size"),blankTheme )
  img = SVG("image/GenDrift_mutation/$N mutation_rate_$μ .svg", 8inch, 6inch)
  draw(img, plt)
  return ev_matrix
end

get_data(gen_drift_p,ev_data,1000)
plot(ev_data, x="x2", Geom.density, Guide.title("Boxplot of mean value fixation occurs"), Guide.XLabel("population size"), blankTheme)
f=convert(Array{Integer}, generation_vector)
plot(ev_matrix, x=ev_matrix[:,485], Geom.histogram, blankTheme)
