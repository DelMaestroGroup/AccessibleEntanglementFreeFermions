# Renyi entanglement entropy of Free fermions in 1D.

push!(LOAD_PATH, joinpath(dirname(@__FILE__), "src"))
using FreeFermionsOEE

using ArgParse
using Printf
using LinearAlgebra

s = ArgParseSettings()
s.autofix_names = true
@add_arg_table s begin
    "L"
        help = "number of sites"
        arg_type = Int
        required = true
    "N"
        help = "number of particles"
        arg_type = Int
        required = true
    "--out"
        metavar = "FILE"
        help = "path to output file"
        required = true
end

add_arg_group(s, "BH parameters")
@add_arg_table s begin
    "--l-min"
        metavar = "l"
        help = "minimum lA"
        arg_type = Int
        default = 1
    "--l-max"
        metavar = "l"
        help = "maximum lA"
        arg_type = Int
        default = 1
    "--l-step"
        metavar = "l"
        help = "lA step"
        arg_type = Int
        default = 1
    "--l-logstep"
        metavar = "Delta"
        help = " log(lA) step"
        arg_type = Float64
        default = 0.1
    "--l-num"
        metavar = "n"
        help = "number of lA"
        arg_type = Int
        default = 1
    "--l-log"
        help = "use logarithmic scale for l"
        action = :store_true
    "--precision"
        metavar = "precision"
        help = "precision"
        arg_type = Int
        default = 54
    "--alpha"
        metavar = "alpha"
        help = "alpha value"
        arg_type = Float64
        default = 2.0
    "--eigs"
        help = "report the spectrum of the projected correlation matrix."
        action = :store_true
    "--prob"
        help = "report the particle-number probability distribution."
        action = :store_true

end

function main()

   c = parsed_args = parse_args(ARGS, s, as_symbols=true)
   # Number of sites
   L = c[:L]
   # Number of particles
   N = c[:N]
   # Output file
   output = c[:out]
   # minimum lA
   lA_min = c[:l_min]
   #precision
   precision=c[:precision]
   alpha=c[:alpha]
   eigenvalues= c[:eigs]
   probabilities= c[:prob]



   if c[:l_log]
      l_cut=trunc(Int,1.0/(exp(c[:l_logstep])-1.0))+1 
      if lA_min<l_cut
         loglnum=trunc(Int,log(c[:l_max]*1.0/l_cut)/c[:l_logstep])
         l_range = zeros(Int, l_cut-lA_min+1+ loglnum)
         for i=1: l_cut-lA_min+1
            l_range[i]=i+ lA_min-1
         end
         for i=1:loglnum
            l_range[i+l_cut-lA_min+1]= trunc(Int,l_cut*exp(i*c[:l_logstep]))
         end
         lA_max=l_range[l_cut+ loglnum-lA_min]
      else
         loglnum=trunc(Int,log(c[:l_max]*1.0/lA_min)/c[:l_logstep])
         l_range = zeros(Int, loglnum+1)
         for i=1: loglnum+1
            l_range[i]= trunc(Int,lA_min*exp((i-1)*c[:l_logstep]))
         end
         lA_max=l_range[loglnum+1]
      end 
   else
      lA_max=(lA_min+c[:l_step]*(c[:l_num]-1))
      l_range = lA_min:c[:l_step]: lA_max
   end
#_______________________
   open(output, "w") do f
      write(f, "# L=$(L), N=$(N), precision=$(precision), alpha=$(alpha) \n")
      write(f, @sprintf "#%7s%24s%24s%24s%24s%24s%24s%24s%24s%24s%24s\n" "l"    "S1"    "S1_OP"    "S_alpha"    "S_alpha_OP(1)"    "sigma2_n" "H_alpha" "navge_alpha" "S_alpha_OP(3)" "S_alpha_OP(4)" "S_alpha_OP(5)")

   for l in l_range
      # ======================= Correlation matrix =============================
      Corr_Matrix = zeros(Float64, l, l)
      phase_factor= Float64
      even_odd=N % 2 
      phase_factor=2*pi/L 
      for i=1: l
         for j=i+1: l
            Corr_Matrix[i,j]=sin(N*phase_factor/2*(i-j))/sin(phase_factor/2*(i-j))/L
            Corr_Matrix[j,i]=Corr_Matrix[i,j]
         end
         Corr_Matrix[i,i]=N/L
      end
	# ==================================End===================================

        # Create C_A
        nu=eigvals!(Symmetric(Corr_Matrix))
        nu_increasing= zeros(Float64, l)
        nu_increasing[l:-1:1]=nu[1:l]
        navge_alpha=0.0 
        for k=1:l
         if nu[k]>=0 && nu[k]<=1
             navge_alpha+=nu[k]^alpha/(nu[k]^alpha+(1-nu[k])^alpha)
         elseif nu[k]>1
             navge_alpha+=1.0
         end
        end

 # Calculate spatial & operational entanglement entropy
        s1_spatial, s1_operational, salpha_spatial, salpha_operational, sigma2_n, LogSumPNalpha, H4, H5, H6 = operational_entanglement(L, N , l, precision, alpha, nu_increasing, eigenvalues, probabilities)
        write(f, @sprintf "%8s%24.12E%24.12E%24.12E%24.12E%24.12E%24.12E%24.12E%24.12E%24.12E%24.12E\n" "$(l)"  s1_spatial s1_operational salpha_spatial salpha_operational sigma2_n LogSumPNalpha navge_alpha  H4 H5 H6)

        flush(f)
   end
   end
end

main()
