
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
l(k) = ((1-θ)*A/2)^(1/(θ+1))*k^(θ/θ+1)

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
function newton(f,initial,ϵ,quietly)
    x = initial
    iter = 1
    while abs(f(x))>= ϵ
        #println("i=$iter val=$x")
        dfdx = ForwardDiff.derivative(f,x)
        if !quietly
            println("dfdx $dfdx")
        end
        x = x - f(x)/dfdx
        iter += 1
    end
    return x
end

#this is a weird answer, I'm not sure what to do about this
@time k1_newton = newton(f,k0,ϵ)
@show k1_newton

# 1.d)
using Plots

#generating guesses for capital, labor, linear from initial to steady state
T = 200 #last time period
ITER = 50 #number of iterations
#guesses for labor and captial
kguesses = [i for i in k0:(kss-k0)/T:kss]
lguesses = [i for i in l0:(lss-l0)/T:lss]

#time series for capital and labor that we're updating
kvec = [i for i in k0:(0.25*kss)/T:kss]
lvec = [i for i in l0:(lss-l0)/T:lss]

#time periods we want to revise
for iter in 1:ITER
    println("iter $iter")
    for t in 2:T
        #FOC for this time period
        f(k) = Euler(c(lvec[t-1],kvec[t-1],k), #consumption of previous period
            lvec[t-1], #labor in previous period
            c(l(k),k,kvec[t+1]), #c_t as a function of k_t
            l(k), #labor in time t as a function of k_t
            kvec[t+1]) #capital in time t+1, a parameter
        kvec[t] = newton(f,kvec[t],ϵ,true)

        if kvec[t] == NaN
            println("leaving")
            break
        end
        lvec[t] = l(kvec[t])
    end
    if iter == 10
        k10 = kvec
    end
end
plot(0:T,kguesses)
plot!(0:T,kvec)

plot(0:T,lguesses)
plot!(0:T,lvec)
