using Gadfly, Distributions, ColorBrewer, DataFrames


# 2-D Evolution Matrix
generation = 500 #  number of generations
repNum= 100      # number of repetitions
genVec = zeros(repNum,1) #timevector t=(x1,x2,...,xN)
EM= zeros(repNum,generation) #evolution matrix

#create dataframe for statistics
evoDf = DataFrame(Any[popSize=Float64[], meanGen=Float64[], std_gen=Float64[], var_gen=Float64[], mute_rate=Float64[]])

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
  println(typeof(M), typeof(N), typeof(μ))


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
            genVec[rep]= gen
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
            genVec[rep]= gen

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

    return gen_drift(), M, genVec
  end

  @time gen_drift(EM, 100, 0.00002)

  #run gen_drift several times for statistics
  ##use try
  ##        catch
  ##function!!!

  function get_data(M::Matrix)      #Float64(x, [, mode::RoundingMode])) FloatRange?

    if nrow(evoDf) >0
      deleterows!(evoDf, 1:nrow(evoDf))
    end
    for N in 100:100:1000
      for μ in 0:10.0^(-10):10.0^(-9)
        gen_drift(EM,N,μ)
        push!(evoDf, [N  μ floor(mean(genVec)) std(genVec) var(genVec)])
      end
    end
println(evoDf)

    plt_data=plot(evoDf, x=evoDf[2], y=evoDf[1], ymin=evoDf[2]-evoDf[3], ymax=evoDf[2]+evoDf[3],
    Geom.point, Geom.line,  Geom.errorbar, Guide.title("meanvalue of generation for extinction/fixation depending on population size"), Guide.XLabel("population size"), Guide.YLabel("generation number"), blankTheme)
    img= PDF("image/GenDrift_mutation/mean.pdf", 8inch, 6inch)
    draw(img, plt_data)
    return draw(img, plt_data)
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

@time simulation(EM, 100, 0.0012)
