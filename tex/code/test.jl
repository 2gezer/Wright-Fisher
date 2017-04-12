using Gadfly
using Distributions
using ColorBrewer
using DataFrames


generation = 500 #  number of generations
rep_nb= 100      # number of repetitions
generation_vector = zeros(rep_nb,1) #timevector t=(x1,x2,...,xN)
EM= zeros(rep_nb,generation) #evolution matrix for genetic drift

EvolutionMatrix= zeros(rep_nb, generation)
@time EM= copy(EvolutionMatrix)
#create dataframe for statistics
ev_data = DataFrame(Any[Pop_size=Float64[], mean_gen=Float64[], std_gen=Float64[], var_gen=Float64[], mute_rate=Float64[]])

#plotting table
blankTheme = Theme(
grid_color=colorant"grey",
panel_fill=colorant"white",
major_label_color=colorant"black",
minor_label_color=colorant"black",
lowlight_opacity=0.9 )


# #simulate pure genetic drift
# function gen_drift_p(M::Matrix, N::Integer)
#
# end

function simulation(M::Matrix, N::Integer, μ::Float64)

  function gen_drift(M::Matrix, N::Integer, μ::Float64)
    row=size(M,1)
    col=size(M,2)
    #mutationrate
    if μ==0
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
      lay =[layer(x=1:generation, y=M[i,:], Geom.path,Theme(default_color=colors[i])) for i in 1:10] #Guide.manual_color_key("legend for plot", ["1","2"]["deepskyblue", "orange"])

      #plot
      plt=plot(lay... , Guide.title("genetic drift"),  Guide.XLabel("generation number"),Guide.YLabel("population number"),blankTheme )
      img = SVG("image/GenDrift_pure/population_size_$N.svg", 8inch, 6inch)
      draw(img, plt)
      return M #compute pure genetic drift

    else
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
      end #genetic drift with mutation

      #plotting window (module @layers using ColorBrewer... color.jl???)
      colors = palette("Spectral", 11)
      lay =[layer(x=1:col, y=M[i,:], Geom.path,Theme(default_color=colors[i])) for i in 1:10] #Guide.manual_color_key("legend for plot", ["1","2"]["deepskyblue", "orange"])
      #plot
      plt=plot(lay... , Guide.title("genetic drift with mutationrate μ = $μ "),  Guide.XLabel("generation number"),Guide.YLabel("population size"),blankTheme )

      #plot histogramm
      plt_h=plot(M,layer(x= M[:,col], Geom.histogram(density=true)),
      layer([x->((N)/(x*(N-x)))],0, N, Geom.line, Theme( default_color=colorant"orange", line_width=2pt)), blankTheme,
      Guide.xlabel("number of A_type individuals (bin size = 10)"),
      Guide.YLabel("density"),
      Guide.title("Histogram of population in generation $col"),
      Guide.manual_color_key("Legend", ["histogram", "theoretical curve"], ["deepskyblue", "orange"]))

      img = PDF("image/GenDrift_mutation/$N mutation_rate_$μ .pdf", 8inch, 6inch) #name of image
      draw(img, vstack(plt,plt_h)) #save both plots in one image
      return M
    end

    return gen_drift(), M, generation_vector
    return @time gen_drift()
  end

 @time gen_drift(EM, 100, 0.0012)

  #run gen_drift several times for statistics


  function get_data(M::Array{Float64,2},  μ::Float64(x, [RoundingMode])) FloatRange
    T=M
    if nrow(ev_data) >0
      deleterows!(ev_data, 1:nrow(ev_data))
    end
    for i in 100:100:1000
      for j in 0:0.00002:0.001
      gen_drift(EM,i,j)
      push!(ev_data, [i  μ floor(mean(generation_vector)) std(generation_vector) var(generation_vector)])
    end
    end

    plt_data=plot(ev_data, x=ev_data[1], y=ev_data[2], ymin=ev_data[2]-ev_data[3], ymax=ev_data[2]+ev_data[3],
    Geom.point, Geom.line,  Geom.errorbar, Guide.title("meanvalue of generation for extinction/fixation depending on population size"), Guide.XLabel("population size"), Guide.YLabel("generation number"), blankTheme)
    if μ==0
      img= PDF("image/GenDrift_pure/mean.pdf", 8inch, 6inch)
    else img= PDF("image/GenDrift_mutation/mean.pdf", 8inch, 6inch)
    end
    draw(img, plt_data)
    return @time get_data()
  end
  @time get_data(EM)
  ##use try
  ##        catch
  ##function!!!

end

# #evaluate data for pure genetic drift
# function get_data_p(df::DataFrame)
#  if nrow(df) >0
#    deleterows!(df, 1:nrow(df))
#  end
#  for i in 100:100:1000
#     gen_drift_m(EM,i,10)
#      push!(df, [i floor(mean(generation_vector)) std(generation_vector)])
#    end
#
#    plt_m=plot(ev_data, x=ev_data[1], y=ev_data[2], ymin=ev_data[2]-ev_data[3], ymax=ev_data[2]+ev_data[3],
#         Geom.point, Geom.errorbar, Guide.title("meanvalue of generation for extinction/fixation depending on population size"), Guide.XLabel("population size"), Guide.YLabel("generation number"), blankTheme)
#    img= SVG("image/GenDrift_pure/mean.svg", 8inch, 6inch)
#    draw(img, plt_m)
#    return gen_drift_m(...);
# end
#return get_data_p()
#for genetic drift with mutation#
#
#     return gen_drift(...), println("image for mutationrate $μ saved")
# end



#simulation(EM,1030,0.042)

@time simulation(EM, 100, 0.0012)
