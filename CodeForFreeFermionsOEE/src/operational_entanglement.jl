"""
Calculate both the spatial and the operational entanglement entropies of a
region A, using the SVD. The latter is the "entanglement of particles"
introduced by Wiseman and Vaccaro in 2003.
"""
function operational_entanglement(L::Int, N::Int, l::Int, precision::Int, alpha::Float64, nu::Vector{Float64}, eigenvalues::Bool, probabilities::Bool)

    setprecision(precision)
#set_bigfloat_precision(precision)
# =======================================================
#Calculating the probability distribution (PN) of the number of particles in the targeted partition.
    scaling_factor= BigFloat
    scaling_factor=10.0^300
    normsDim=l+1
    #P1 = zeros(BigFloat, 2)
    PN_old = zeros(BigFloat, normsDim)
    PN_new = zeros(BigFloat, normsDim)
    if nu[1]>=0 && nu[1]<=1
       nuVal=nu[1]
    elseif nu[1]>1
       nuVal=1.0
    else
       nuVal=0.0
    end
    PN_old[1]=scaling_factor*(1-nuVal)                         
    PN_old[2]=scaling_factor*nuVal
    NN=2
    for i= 2:normsDim-1
       if nu[i]>=0 && nu[i]<=1
          nuVal=nu[i]
       elseif nu[i]>1
          nuVal=1.0
       else
          nuVal=0.0
       end
       PN_new[1]=(1-nuVal)*PN_old[1] 
       PN_new[NN+1]=nuVal*PN_old[NN] 
       for k= 1:NN-1
          PN_new[k+1]=(1-nuVal)*PN_old[k+1]+nuVal*PN_old[k]
       end
       NN=NN+1
       for k=1: NN
          PN_old[k]=PN_new[k]
       end
    end
    if l==1
       for k=1: normsDim
          PN_new[k]=PN_old[k]
       end
    end
    PN_new/= scaling_factor
    err_PN = abs(sum(PN_new) - 1.0)
    if err_PN > 1e-12
        @warn("probability distribution error: $(err_PN)")

    end
# =======================================================
# =======================================================
#Calculating the second Renyi entanglement entropy.
    Salpha = zeros(BigFloat, 3)
    Salpha_old = zeros(BigFloat, normsDim)
    Salpha_new = zeros(BigFloat, normsDim)
    #alpha=2.0
    if nu[1]>=0 && nu[1]<=1
       nuVal=nu[1]
    elseif nu[1]>1
       nuVal=1.0
    else
       nuVal=0.0
    end
    Salpha_old[1]=scaling_factor*(1-nuVal)^alpha
    Salpha_old[2]=scaling_factor*nuVal^alpha
    NN=2
    for i= 2:normsDim-1
       if nu[i]>=0 && nu[i]<=1
          nuVal=nu[i]
       elseif nu[i]>1
          nuVal=1.0
       else
          nuVal=0.0
       end
       Salpha_new[1]=(1-nuVal)^alpha*Salpha_old[1]
       Salpha_new[NN+1]=nuVal^alpha*Salpha_old[NN]
       for k= 1:NN-1
          Salpha_new[k+1]=(1-nuVal)^alpha*Salpha_old[k+1]+nuVal^alpha*Salpha_old[k]
       end
       NN=NN+1
       for k=1: NN
          Salpha_old[k]=Salpha_new[k]
       end
    end


    if l==1
       for k=1: normsDim
          Salpha_new[k]=Salpha_old[k]
       end
    end
    Salpha_new/= scaling_factor
# =======================================================
    LogSumPNalpha= BigFloat
    SalphalogSalpha= BigFloat
    PnlogPN= BigFloat
    sigma2_n= BigFloat
    sigma2_n_PN= BigFloat
    n= BigFloat
    n2= BigFloat

    H4= BigFloat
    H5= BigFloat
    H6= BigFloat

    S1_sp= BigFloat 
    Salpha_sp= BigFloat 
    S1_op= BigFloat 
    Salpha_op= BigFloat
# =======================================================

    S1_sp=0.0
    Salpha_sp=0.0
    for i= 1:normsDim-1
       if nu[i]>0 && nu[i]<1
             S1_sp-=(nu[i]*log(nu[i])+(1-nu[i])*log((1-nu[i])))
             Salpha_sp-=log(nu[i]^alpha+(1-nu[i])^alpha) 
       end
    end
    err_Salpha = abs(sum(Salpha_new) - exp(-1.0*Salpha_sp))/exp(-1.0*Salpha_sp)
    if err_Salpha > 1e-12
        @warn("Salpha_op error: $(err_Salpha)")

    end
    Salpha_sp/=(alpha-1)
# =======================================================

    sigma2_n=0.0
    for i= 1:normsDim-1
       if nu[i]>0 && nu[i]<1
          sigma2_n+=(nu[i]*(1-nu[i]))
       end
    end

    sigma2_n_PN=0.0
    n=0
    n2=0
    for i= 1: normsDim
       if PN_new[i]>0
          n+=PN_new[i]*(i-1)
          n2+=PN_new[i]*(i-1)^2
       end
    end
    sigma2_n_PN=n2-n^2

    err_sigma2 = abs(sigma2_n_PN - sigma2_n)/sigma2_n
    if err_sigma2 > 1e-12
        @warn("sigma2_n error: $(err_sigma2)")

    end
# =======================================================

    PnlogPN=0.0
    for i= 1: normsDim
       if PN_new[i]>0
          PnlogPN+=PN_new[i]*log(PN_new[i])
       end
    end
# =======================================================



    S1_op=0.0
    S1_op=S1_sp+PnlogPN

    PNlogSalpha=0.0
    for i= 1: normsDim
       if Salpha_new[i]>0
          PNlogSalpha-=PN_new[i]*log(Salpha_new[i])
       end
    end

    Salpha_op= PNlogSalpha+alpha*PnlogPN
    Salpha_op/=(alpha-1)

    LogSumPNalpha=sum(PN_new.^alpha)
    LogSumPNalpha=-log(LogSumPNalpha)/(alpha-1)

# =======================================================

    H4=0.0
    H5=0.0
    H6=0.0
    H = zeros(BigFloat, normsDim)
    for i= 1: normsDim
       if PN_new[i]>0
          H[i]+=Salpha_new[i]/PN_new[i]^alpha
          H4+=PN_new[i]*H[i]
          H5+=PN_new[i]*H[i]^(1/(alpha-1))
          H6+=Salpha_new[i]^(1/alpha)
       end

    end

    H4=log(H4)/(1-alpha)
    H5=-log(H5)
    H6=-log(H6)*(alpha/(alpha-1))

# =======================New=============================
# =======================================================



     if l==100
        if probabilities
           Output="L$(L)N$(N)l$(l)alpha$(alpha)Pn.dat"
           open(Output, "w") do fp
              write(fp, "# L=$(L), N=$(N), l=$(l), alpha=$(alpha) \n")
              write(fp, "# ==============================================================================\n") 
              write(fp, @sprintf "#%7s%24s%24s%24s\n" "n"  "P(n)" "P(n)^alpha"  "P(n)_alpha" )

              sumPnalpha=0.0
              sumPn_alpha=0.0
              sumSalpha_new =sum(Salpha_new)
              sumPN_newalpha =sum(PN_new.^alpha)
              for i=1: normsDim
                 write(fp, @sprintf "%8s%24.12E%24.12E%24.12E\n" "$(i-1)"  PN_new[i]  PN_new[i]^alpha/sumPN_newalpha Salpha_new[i]/sumSalpha_new)
                 sumPnalpha= sumPnalpha+PN_new[i]^alpha/sumPN_newalpha
                 sumPn_alpha= sumPn_alpha+Salpha_new[i]/sumSalpha_new
              end
               write(fp, "#$(sumPnalpha)  \n")
               write(fp, "#$(sumPn_alpha)  \n")
               flush(fp)
           end
        end

#___________________________________________eigvals_______________________________________________
        if eigenvalues
           Outputn="L$(L)N$(N)l$(l)alpha$(alpha)eigvals.dat"
           open(Outputn, "w") do fn
              write(fn, "# L=$(L), N=$(N), l=$(l) \n")
              write(fn, "# ===============================\n") 
              write(fn, @sprintf "#%7s%24s\n" "n"  "nu" )
              for i=1: normsDim-1
                 if nu[i]>=0 && nu[i]<=1
                    nuVal=nu[i]
                 elseif nu[i]>1
                    nuVal=1.0
                 else
                    nuVal=0.0
                 end
                 write(fn, @sprintf "%8s%24.12E\n" "$(i)"  nuVal)
              end
              flush(fn)
           end
        end
#__________________________________________________________________________________________

     end


    S1_sp, S1_op, Salpha_sp, Salpha_op, sigma2_n,LogSumPNalpha,H4,H5,H6
end
