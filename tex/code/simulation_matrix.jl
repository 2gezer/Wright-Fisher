# packages
  using Gadfly
  using Distributions
  using DataFrames


#work with methods


  # size of future matrix
  generation::Integer = 500 # ==x
  rep_nb::Integer= 100  # ==y  number of repetitions

 #timevector t=(x1,x2,...,xN)
  generation_vector = zeros(rep_nb,1)
  #evolution matrix
  ev_matrix= zeros(rep_nb,generation)


#simulate pure genetic drift
function gen_drift(N,number)

#

  f(N)= ev_matrix[:,1] = N/2
   for rep in 1:repetition
     for gen in 2:generation
        #next generation will have random number of individuals A with probability of current generation
        p= ev_matrix[rep,gen-1]/N
        A= rand(Binomial(N,p))

        #update state ⃗
        ev_matrix[rep,gen]=A

        #save the generation number and break if extinction/fixation occurs
        if ev_matrix[rep,gen]==0 || ev_matrix[rep,gen]==N
          generation_vector[rep]= gen
          ev_matrix[rep,gen+1:generation]= NaN
          break
        end
      end
    end
  end
  #plot
for i in 1: number
ga
layers = [layer(x=t, y=x[ii,:],Geom.path,Theme(default_color=colors[ii]))
    for ii in 1:11];
  # m_gen =[]
  # err_gen=[]
  # pop_size=collect(100:100:1000)

# #plot population evolution for 5 repetitions
# for i in 50:55
# plot(ev_matrix[i,:]/100)
# title("genetic drfit with mutation rate 0.001")
# xlabel("generation number")
# ylabel("A-type individuals")
# end
plt=
#plot meanvalue of the generations, when idiviual takes over, depending on populationsize
# #plot mean value extinction/fixation against population size
# errorbar(pop_size,m_gen,err_gen)
# title("mean value of extinction/fixation dependending on Population size")
# xlabel("Population size N")
# ylabel("generation")

# #plot histogram of repetition columns
# min_t= minimum(floor(Int64,generation_vector))
# max_t= maximum(floor(Int64,generation_vector))
# h = PyPlot.plt[:hist](ev_matrix[:,max_t])
# title(["Histogram of evolution of A-type in generation: " max_t])
# xlabel("population number")
# ylabel("A-type individuals ")

end

function simulation(N)
  ev_matrix[:,1] = N/2
   for rep in 1:repetition
     for gen in 2:generation

       #add mutation
       μ=0.001
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
  end

simulation(100)
