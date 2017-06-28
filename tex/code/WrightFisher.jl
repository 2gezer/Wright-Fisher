using Gadfly, Distributions, ColorBrewer, DataFrames


# 2-D Evolution Matrix
generation = 500 #  number of generations
<<<<<<< HEAD
rep_nb= 100      # number of repetitions
generation_vector = zeros(rep_nb,1) #timevector t=(x1,x2,...,xN)
EM= zeros(rep_nb,generation) #evolution matrix for genetic drift
=======
repNum= 100      # number of repetitions
genVec = zeros(repNum,1) #timevector t=(x1,x2,...,xN)
EM= zeros(repNum,generation) #evolution matrix
>>>>>>> 3e7f40d6b5aed01069dd14df00e572b0b929d6d4

EvolutionMatrix= zeros(rep_nb, generation)
@time EM= copy(EvolutionMatrix)
#create dataframe for statistics
<<<<<<< HEAD
ev_data = DataFrame(Any[Pop_size=Float64[], mean_gen=Float64[], std_gen=Float64[], var_gen=Float64[], mute_rate=Float64[]])
=======
evoDf = DataFrame(Any[popSize=Float64[], A_types=Integer[], meanGen=Float64[], std_gen=Float64[], var_gen=Float64[], mute_rate=Float64[]])
>>>>>>> 3e7f40d6b5aed01069dd14df00e572b0b929d6d4

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
<<<<<<< HEAD
=======
  println(typeof(M), typeof(N), typeof(μ))

>>>>>>> 3e7f40d6b5aed01069dd14df00e572b0b929d6d4

  function gen_drift(M::Matrix, N::Integer, μ::Float64)
    row=size(M,1)
    col=size(M,2)
    #mutationrate
    if μ==0
      M[:,1] = N/2
      for rep in 1:row
        for gen in 2:col
          #next generation will have random number of individuals A with probability of current generation
          p= M[rep,gen-1]/N
          A= rand(Binomial(N,p))

          #update state
          M[rep,gen]=A

          #save the generation number and break if extinction/fixation occurs
          if M[rep,gen]==0 || M[rep,gen]==N
<<<<<<< HEAD
            generation_vector[rep]= gen
=======
            genVec[rep]= gen
>>>>>>> 3e7f40d6b5aed01069dd14df00e572b0b929d6d4
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
<<<<<<< HEAD
            generation_vector[rep]= gen
=======
            genVec[rep]= gen
>>>>>>> 3e7f40d6b5aed01069dd14df00e572b0b929d6d4

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
<<<<<<< HEAD

    return gen_drift(), M, generation_vector
    return @time gen_drift()
  end

 @time gen_drift(EM, 100, 0.00002)

  #run gen_drift several times for statistics
  ##use try
  ##        catch
  ##function!!!

  function get_data(M::Array{Float64,2}, N::Integer, μ::Float64)      #Float64(x, [, mode::RoundingMode])) FloatRange?

    if nrow(ev_data) >0
      deleterows!(ev_data, 1:nrow(ev_data))
    end
    for N in 100:100:1000
      for μ in 0:10.0^(-10):10.0^(-9)
      gen_drift(EM,N,μ)
      push!(ev_data, [N  μ floor(mean(generation_vector)) std(generation_vector) var(generation_vector)])
    end
    end

    plt_data=plot(ev_data, x=ev_data[1], y=ev_data[2], ymin=ev_data[2]-ev_data[3], ymax=ev_data[2]+ev_data[3],
    Geom.point, Geom.line,  Geom.errorbar, Guide.title("meanvalue of generation for extinction/fixation depending on population size"), Guide.XLabel("population size"), Guide.YLabel("generation number"), blankTheme)
    if μ==0
      img= PDF("image/GenDrift_pure/mean.pdf", 8inch, 6inch)
    else img= PDF("image/GenDrift_mutation/mean.pdf", 8inch, 6inch)
    end
    draw(img, plt_data)
    return get_data()
  end

get_data


  try
    @time get_data(EM)
    println("Test")
catch ErrorException
end
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

=======

    return gen_drift(), M, genVec
  end

  @time gen_drift(EM, 100, 0.00002)

  #run gen_drift several times for statistics
  ##use try
  ##        catch
  ##function!!!

  function get_data(M::Matrix)      #Float64(x, [, mode::RoundingMode])) FloatRange?
    c=Integer(floor(mean(genVec)))
    d=Integer(floor(mean(EM[:,c])))
    if nrow(evoDf) >0
      deleterows!(evoDf, 1:nrow(evoDf))
    end
    for N in 100:100:1000
      for μ in 10.0^(-7):10.0^(-10):10.0^(-9)
        gen_drift(EM,N,μ)
        push!(evoDf, [N d μ c std(genVec) var(genVec)])
        println(μ)
      end
    end

writetable("output.txt", evoDf, separator = ',', header = true)

    plt_data=plot(evoDf, x=evoDf[1], y=evoDf[1], ymin=evoDf[2]-evoDf[3], ymax=evoDf[2]+evoDf[3],
    Geom.point, Geom.line,  Geom.errorbar, Guide.title("meanvalue of generation for extinction/fixation depending on population size"), Guide.XLabel("population size"), Guide.YLabel("generation number"), blankTheme)
    img= PDF("image/GenDrift_mutation/mean.pdf", 8inch, 6inch)
    draw(img, plt_data)
    return get_data()
  end

get_data(EM)
  println("result")

  # try
  #   @time get_data(EM)
  #
  # catch ErrorException
  #   println("get_data throughs ErrorException :(")
  # end

end

>>>>>>> 3e7f40d6b5aed01069dd14df00e572b0b929d6d4
@time simulation(EM, 100, 0.0012)
