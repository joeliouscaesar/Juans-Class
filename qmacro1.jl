
#parameters
σ = 2
β = 0.99
θ = 0.3
δ = 0.07
A = 0.6

#1.a)

#we did all the math this just calculates kss and lss

#this is k_ss/l_ss
klratio = ((1-β*(1-δ))/(β*θ*A))^(1/(θ-1))

lss = (1-θ)*A/2*(klratio)^θ
kss = lss*klratio

#1.b)

#bisection solver!
#f is function
#left and right are domain endpoints
#ϵ is an error tolerance
#setting quietly to true mutes output
function bisection(f,left,right,ϵ,quietly=false)
    if f(left)*f(right) > 0
        println("Invalid domain endpoints")
        return NaN
    end

    middle = (left+right)/2
    iter = 1
    while abs(f(middle)) >= ϵ
        if !quietly
            println("i=$iter val=$middle")
        end

        if sign(f(middle))== sign(f(left))
            left = middle
        else
            right = middle
        end
        middle = (left+right)/2
        iter+=1
    end
    return middle
end

#writing labor as a function of capital, from labor FOC
l(k) = ((1-θ)*A/2)^(1/(θ+1))*k^(θ/(θ+1))

#Euler set equal to zero, we'll apply bisection to this
#I did RHS to make it more readable?
Euler_RHS(ct1,lt1,kt1) = β*(ct1-lt1^2)^(-σ)*(θ*A*(lt1/kt1)^(1-θ)+1-δ)
Euler(ct,lt,ct1,lt1,kt1) = (ct-lt^2)^(-σ) - Euler_RHS(ct1,lt1,kt1)

#getting l0
k0 = 0.75*kss
l0 = l(k0)
k2 = 0.85*kss

@show k0
@show l0
@show k2

#c as a function of lt, kt and kt1 from constraint
c(lt,kt,kt1) = A*(kt)^(θ)*(lt)^(1-θ)+(kt)*(1-δ)-(kt1)

#function we'll be doing bisection on
f(k1) = Euler(c(l0,k0,k1), #this is consumption in time zero, function of k1
    l0, #labor in time 0, solved from k0
    c(l(k1),k1,k2), #consumption in time 1 as a function of k1
    l(k1), #labor in time 1 as a function of k1
    k2) #capital in time 2, a parameter

#thought k0 and k2 made sense as endpoints, it's really sensitive
@show f(k0)
@show f(k2)
#error tolerance
ϵ = 10e-7
@show ϵ
@time k1 = bisection(f,k0,k2,ϵ)
@show k1

# 1.c)

using ForwardDiff
#newton solver

#ϵ is the error tolerance
function newton(f,initial,ϵ,quietly=false)
    x = initial
    iter = 1
    while abs(f(x))>= ϵ
        #println("i=$iter val=$x")
        dfdx = ForwardDiff.derivative(f,x)
        if !quietly
            println("i=$iter val=$x")
        end
        x = x - f(x)/dfdx
        iter += 1
    end
    return x
end

#this is a weird answer, I'm not sure what to do about this
@show k0
@show k2
@show mid = (k0+k2)/2
@time k1_newton = newton(f,(k0+k2)/2,ϵ)
@show k1_newton


# 1.d)
using Plots

#function which calculates the captial accumulation path using Newton's method
#T is the last time period
#ITER is the number of iterations
#kss is the steady state capital, k0 is the initial value
#lss is the steady state labor, l0 is the initial value
function capital_path(T,ITER,kss,k0,lss,l0,quietly=false)

    #guesses for labor and capital, straight line to steady state
    kguesses = [i for i in k0:(kss-k0)/T:kss]
    lguesses = [i for i in l0:(lss-l0)/T:lss]

    #time series for capital and labor that we're updating
    kvec = copy(kguesses)
    lvec = copy(lguesses)

    for iter in 1:ITER
        if !quietly
            println("iter $iter")
        end
        for t in 2:T
            #gets the FOC for this period
            f(k) = Euler(c(lvec[t-1],kvec[t-1],k), #consumption of previous period
                lvec[t-1], #labor in previous period
                c(l(k),k,kvec[t+1]), #c_t as a function of k_t
                l(k), #labor in time t as a function of k_t
                kvec[t+1]) #capital in time t+1, a parameter
            #updating capital and labor
            kvec[t] = newton(f,kvec[t],ϵ,true)
            lvec[t] = l(kvec[t])
        end
    end
    return kvec
end

#plots the initial guess line and after 50 iterations
#for capital and labor
T = 50

#capital paths for different numbers of iterations
#j0 is zero iterations, j10 is 10
kvec_j0 = capital_path(T,0,kss,k0,lss,l0)
kvec_j1 = capital_path(T,1,kss,k0,lss,l0)
kvec_j2 = capital_path(T,2,kss,k0,lss,l0)
kvec_j5 = capital_path(T,5,kss,k0,lss,l0)
kvec_j10 = capital_path(T,10,kss,k0,lss,l0)
kvec_j50 = capital_path(T,50,kss,k0,lss,l0)
kvec_j5000 = capital_path(T,5000,kss,k0,lss,l0)
#50000 looks same as 5000, takes about 2 minutes
#@time kvec_j50000 = capital_path(T,50000,kss,k0,lss,l0,true)

#capital
plot(0:T,kvec_j0,label="j=0",legend = :outertopleft,title="Capital Path")
plot!(0:T,kvec_j1,label="j=1")
plot!(0:T,kvec_j2,label="j=2")
plot!(0:T,kvec_j5,label="j=3")
plot!(0:T,kvec_j10,label="j=4")
plot!(0:T,kvec_j50,label="j=5")
plot!(0:T,kvec_j5000,label="Optimal")

# 1.5)
invest = kvec_j5000
labor = l.(kvec_j5000)

output(k,l) = A*k^θ*l^(1-θ)+k*(1-δ)
outputpath = [output(invest[t],labor[t]) for t in 1:51]

#this is only 50...
consumption = [c(labor[t],invest[t],invest[t+1]) for t in 1:50]

#adds steady state consumption to the consumption path
push!(consumption,c(lss,kss,kss))

plot(0:T,invest,label="investment",legend = :outertopleft)
plot!(0:T,consumption,label="consumption")
plot!(0:T,labor,label="employment")
plot!(0:T,outputpath,label="output")

#quick consistency check that next period investment plus current consumption
#is equal to output
invcons = [invest[t+1]+consumption[t] for t in 1:50]
push!(invcons,c(lss,kss,kss)+kss)
invcons-outputpath
@show sum(invcons-outputpath)
