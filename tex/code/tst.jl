using Gadfly
using Distributions
using ColorBrewer
using DataFrames


generation = 500 #  number of generations
rep_nb= 100      # number of repetitions
generation_vector = zeros(rep_nb,1) #timevector t=(x1,x2,...,xN)
matrix_pure= zeros(rep_nb,generation) #evolution matrix for genetic drift
matrix_mutation= zeros(100,1000)     #evolution matrix for genetic drift with mutation

#create dataframe for statistics
ev_data = DataFrame(Any[Pop_size=Float64[], m_genv=Float64[], err_gen=Float64[]])

#plotting table
blankTheme = Theme(
grid_color=colorant"grey",
panel_fill=colorant"white",
major_label_color=colorant"black",
minor_label_color=colorant"black")


#simulate pure genetic drift
function gen_drift_p(M::Matrix, N::Integer, number::Integer)
  M[:,1] = N/2
  for rep in 1:rep_nb
    for gen in 2:generation
      #next generation will have random number of individuals A with probability of current generation
      p= M[rep,gen-1]/N
      A= rand(Binomial(N,p))

      #update state
      M[rep,gen]=A

      #save the generation number and break if extinction/fixation occurs
      if M[rep,gen]==0 || M[rep,gen]==N
        generation_vector[rep]= gen
        M[rep,gen+1:generation]= NaN
        break
      end
    end
  end
  #plotting window (module @layers using ColorBrewer using Gadfly end...???)
  colors = palette("Spectral", 11)
  lay =[layer(x=1:generation, y=M[i,:], Geom.path,Theme(default_color=colors[i])) for i in 1:number] #Guide.manual_color_key("legend for plot", ["1","2"]["deepskyblue", "orange"])

  #plot
  plt=plot(lay... , Guide.title("genetic drift"),  Guide.XLabel("generation number"),Guide.YLabel("population number"),blankTheme )
  img = SVG("image/GenDrift_pure/population_size_$N.svg", 8inch, 6inch)
  draw(img, plt)
  return M
end

function gen_drift_m(M::Matrix, N::Integer,number::Integer)
  row=size(M,1)
  col=size(M,2)
  μ=0.001 #mutationrate
  M[:,1] = N/2
  for rep in 1:row
    for gen in 2:col

      #add mutation
      a=M[rep,gen-1]
      a_p=rand(Poisson(μ*a))
      b=N-M[rep,gen-1]
      b_p=rand(Poisson(μ*b))
      mutate = b_p-a_p
      M[rep,gen]=M[rep,gen-1]+mutate

      #genetic drift loop
      p= M[rep,gen]/N
      A= rand(Binomial(N,p))

      #update state
      M[rep,gen]=A

      #count extinction/fixation
      if M[rep,gen]==0 || M[rep,gen]==N
        generation_vector[rep]= gen

      end
    end
  end

  #plotting window (module @layers using ColorBrewer... color.jl???)
  colors = palette("Spectral", 11)
  lay =[layer(x=1:col, y=M[i,:], Geom.path,Theme(default_color=colors[i])) for i in 1:number] #Guide.manual_color_key("legend for plot", ["1","2"]["deepskyblue", "orange"])
  #plot
  plt=plot(lay... , Guide.title("genetic drift with mutationrate μ = $μ "),  Guide.XLabel("generation number"),Guide.YLabel("population size"),blankTheme )

  #plot histogramm
  plt_h=plot(M,layer(x= M[:,col], Geom.histogram(bincount=50,density=true)),
             layer([x->(N/(x*(N-x)))], 0, N, Geom.line, Theme(default_color=colorant"orange", key_title_color=colorant"orange")), blankTheme,
             Guide.xlabel("population size"),
             Guide.title("Histogram of population in generation $col"),
             Guide.manual_color_key("Legend", ["histogram", "theoretical curve"], ["deepskyblue", "orange"]))

  img = SVG("image/GenDrift_mutation/$N mutation_rate_$μ .svg", 8inch, 6inch) #name of image
  draw(img, vstack(plt,plt_h)) #save both plots in one image

  return M
end

#evaluate data for genetic drift with mutation
function get_data_p(df::DataFrame)
 if nrow(df) >0
   deleterows!(df, 1:nrow(df))
 end
 for i in 100:100:1000
    gen_drift_m(matrix_mutation,i,10)
     push!(df, [i floor(mean(generation_vector)) std(generation_vector)])
   end

   plt_m=plot(ev_data, x=ev_data[1], y=ev_data[2], ymin=ev_data[2]-ev_data[3], ymax=ev_data[2]+ev_data[3],
        Geom.point, Geom.errorbar, Guide.title("meanvalue of generation for extinction/fixation depending on population size"), Guide.XLabel("population size"), Guide.YLabel("generation number"), blankTheme)
   img= SVG("image/Gendrift_mutation/mean.svg", 8inch, 6inch)
   draw(img, plt_m)
   return df
end

#for pure gentic drift
function get_data_m(df::DataFrame)
 if nrow(df) >0
   deleterows!(df, 1:nrow(df))
 end
 for i in 100:100:1000
    gen_drift_p(matrix_pure,i,10)
     push!(df, [i floor(mean(generation_vector)) std(generation_vector)])
   end

   plt_m=plot(ev_data, x=ev_data[1], y=ev_data[2], ymin=ev_data[2]-ev_data[3], ymax=ev_data[2]+ev_data[3],
        Geom.point, Geom.errorbar, Guide.title("meanvalue of generation for extinction/fixation depending on population size"), Guide.XLabel("population size"), Guide.YLabel("generation number"), blankTheme)
   img= SVG("image/Gendrift_pure/mean.svg", 8inch, 6inch)
   draw(img, plt_m)
   return df
end

get_data(ev_data)

plot(ev_data, x="x2", Geom.density, Guide.title("Boxplot of mean value fixation occurs" , Guide.XLabel("population size"), blankTheme)
f=convert(Array{Integer}, generation_vector)

gen_drift_m(matrix_mutation,100,10)
