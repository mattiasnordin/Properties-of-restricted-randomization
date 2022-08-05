#This code has been tested with Julia version 1.7.2
using Combinatorics #Tested with version 1.0.2
using Distributions #Tested with version 0.23.12
using OrderedCollections #Tested with version 1.4.1

function DM_est(vec)
    #This function takes an assignment vector as input and calculate the
    #difference-in-means estimate (equation 1)
    return (transpose(vec) * Y1 .- transpose(one_vec .- vec) * Y0)[1] / N2
end

#Set sample size (N≥8) and number of assignment vectors in each design (H≥4).
#Note that the code is not generated to be
#efficient but rather transparent and consistent with the paper. The
#computational complexity of the code increases very fast in N and H.
#It is therefore only computationally feasible to work with small values
#of N and H. Both N and H have to be even numbers.
N = 8
N2 = Int(N/2)
H = 4
H2 = Int(H/2)
Nₐ = binomial(N, N2)
Nₐ2 = Int(Nₐ/2)

#Generate all the assignment vectors and store them (lexicographically)
#in the ordered set 𝒲. This is equivalent to the complete randomization design
𝒲 = OrderedSet()
for j = combinations(1:N, N2)
    w = reshape(Int.(zeros(N)), N, 1)
    w[j] .+= 1
    push!(𝒲, w)
end

#Generate some data. The data can be anything and is not restricted to
#any particular distribution
Y0 = reshape(rand(Normal(0, 1), N), N, 1)
Y1 = Y0 .+ rand(Normal(0, 1), N)
τ = mean(Y1) - mean(Y0)

#Generate the Nₐ different treatment effect estimates
one_vec = reshape(Int.(ones(N)), N, 1)
τ_hat = Array{Float64}(undef, Nₐ)
for j = 1:Nₐ
    τ_hat[j] = DM_est(𝒲[j])
end

#Equation 7 gives the MSE of complete randomization
σ² = 1 / Nₐ * sum((τ_hat .- τ).^2)

#Generate the full set of designs of size H, 𝒦ₕ. Go through the first
#half of the lexicographic ordering and add the mirrors from the second half.
𝒦ₕ = Set()
for i = combinations(1:Nₐ2, H2)
    𝒲ₕᵏ = Set()
    for j = i
        push!(𝒲ₕᵏ, 𝒲[j])
        push!(𝒲ₕᵏ, 𝒲[Nₐ+1-j])
    end
    push!(𝒦ₕ, 𝒲ₕᵏ)
end

#Equation 10 gives the variance of the MSE for the set 𝒦ₕ
VarMSE𝒦ₕ = 1/binomial(Nₐ2,H2)*
           sum([(1/H*sum([(DM_est(wⱼ) - τ)^2 for wⱼ = 𝒲ₕᵏ]) -  σ²)^2
                for 𝒲ₕᵏ = 𝒦ₕ])

#Equation A8
m = binomial(Nₐ2, H2)
VarMSE𝒦ₕ_eqA8 = 1/m * sum([
    σ²^2 - 2 / H *
    sum([(DM_est(wⱼ) - τ)^2 for wⱼ = 𝒲ₕᵏ]) * σ² + 1 / H^2 *
    sum([(DM_est(wⱼ) - τ)^2 * (DM_est(wⱼᵖ) - τ)^2 for wⱼ = 𝒲ₕᵏ, wⱼᵖ = 𝒲ₕᵏ])
for 𝒲ₕᵏ = 𝒦ₕ])

#Generate the Nₐ array of squared difference between estimate and τ
r = (τ_hat .- τ).^2
#Generate the useful sums of r
sum_r2 = sum(r.^2)
sum_rr = sum([r[j]*r[jp] for j=1:Nₐ, jp=1:Nₐ if (j!=jp) & (j!=Nₐ+1-jp)])

#Equation A9
VarMSE𝒦ₕ_eqA9 = σ²^2 - 2 * σ² / Nₐ * sum(r) +
                2 / (H * Nₐ) * sum_r2 + (H-2)/(H * Nₐ * (Nₐ-2)) * sum_rr

#Generate the left hand side and right hand side of equation A10.
#Assert that these are equivalent (except for rounding error due
#to floating-point arithmetic).
lhs_eqA10 = σ²^2 - 2 * σ² / Nₐ * sum(r)
rhs_eqA10 = -2 / Nₐ^2 * sum_r2 - 1 / Nₐ^2 * sum_rr
@assert abs(lhs_eqA10-rhs_eqA10) < 1e-12

#Equation A11
p̄ = sum_r2 / Nₐ
q̄ = 1 / (Nₐ * (Nₐ - 2)) * sum_rr
VarMSE𝒦ₕ_eqA11 = 2 * (Nₐ - H) / (H * Nₐ) * (p̄ - q̄)

#𝒦ₕt⊆𝒦ₕ which satisfies conditions 1 and 2. The code below generates
#an example of such a set. This particular example generates
#a set where all assignment vectors in any design in
#𝒦ₕt have either uniqueness u or N/2-u where 1 ≤ u ≤ N/2-1, except for
#all the mirrors who have uniqueness N/2. This set could be changed to
#any other set which satisfies conditions 1 and 2. For instance, 𝒦ₕ
#itself is such a set.
𝒦ₕt = Set()
u = 1
for i = combinations(1:Nₐ2, H2)
    bool_array = [(sum(𝒲[i[j]] .- 𝒲[i[jᵖ]] .== 1) == u) |
                  (sum(𝒲[i[j]] .- 𝒲[i[jᵖ]] .== 1) == N2-u)
                  for j = 1:H2-1 for jᵖ=j+1:H2]
    if minimum(bool_array) == true
        𝒲ₕᵏ = Set()
        for j = i
            push!(𝒲ₕᵏ, 𝒲[j])
            push!(𝒲ₕᵏ, 𝒲[Nₐ+1-j])
        end
        push!(𝒦ₕt, 𝒲ₕᵏ)
    end
end

#Cardinality of 𝒦ₕt
m = length(𝒦ₕt)

#The variance of the MSE from the set 𝒦ₕt is
VarMSE𝒦ₕt = 1/m*sum([(1/H*sum(
    [(DM_est(wⱼ) - τ)^2 for wⱼ = 𝒲ₕᵏt]) -  σ²)^2 for 𝒲ₕᵏt = 𝒦ₕt])

#Create vᵤ𝒦ₕt in line with equation A13
u_list = [[sum(wⱼ-wⱼᵖ.==1) for wⱼ = 𝒲ₕᵏ, wⱼᵖ = 𝒲ₕᵏ] for 𝒲ₕᵏ = 𝒦ₕt]
u_list = [(u_list...)...]
vᵤ𝒦ₕt = Dict([(u, 1 / (m * H^2) * count(x->x==u, u_list)) for u = 0:N2])

#Generate nᵤ
nᵤ = Dict([(u, Nₐ * binomial(N2, u) * binomial(N2, N2-u)) for u = 0:N2])

#Generate equation A14 by first generating a vector of length Nₐ² with
#all pairwise uniquenesses and another with all combinations of ra, rb.
#Then sum all ra, rb with a given u in rr_sum_all. From this, equation A14
#can be generated.
u_all = vec([sum(𝒲[a]-𝒲[b].==1) for a = 1:Nₐ, b = 1:Nₐ])
rr_all = vec([(DM_est(𝒲[a]) - τ)^2 * (DM_est(𝒲[b]) - τ)^2
              for a = 1:Nₐ, b = 1:Nₐ])
rr_sum_all = Dict([(u, sum([rr_all[j] for j = 1:Nₐ^2 if u_all[j] == u]))
                   for u=0:N2])
VarMSE𝒦ₕt_eqA14 = sum([vᵤ𝒦ₕt[u] / nᵤ[u] * rr_sum_all[u] for u = 0:N2]) -
                  1 / Nₐ^2 * sum(rr_all)

#Generate q̄ᵤ
q̄ᵤ = Dict([(u, 1 / nᵤ[u] * rr_sum_all[u]) for u = 0:N2])

#Equation A15
VarMSE𝒦ₕt_eqA15 = sum([(vᵤ𝒦ₕt[u] - nᵤ[u] / Nₐ^2) * q̄ᵤ[u] for u = 0:N2])

α = []
for i = 1:Nₐ
    alpha = deepcopy(Y1)
    alpha[𝒲[i].==0] .= -Y0[𝒲[i].==0]
    push!(α, alpha)
end

#Equation A18
rarb = []
for a = 1:Nₐ
    for b = 1:Nₐ
        push!(rarb,
               τ^4 -
               4*τ^3/N*(
                    sum([α[a][i] + α[b][i] for i = 1:N])
               ) +
               4*τ^2/N^2*(
                    sum([α[a][i]^2 + α[b][i]^2 for i = 1:N]) +
                    4*sum([α[a][i]*α[b][j] for i = 1:N for j = 1:N]) +
                    2*sum([α[a][i]*α[a][j] + α[b][i]*α[b][j]
                         for i = 1:N-1 for j = i+1:N])
               ) -
               16*τ/N^3*(
                    sum([α[a][i]^2*α[b][j] + α[a][i]*α[b][j]^2
                         for i = 1:N for j = 1:N]) +
                    2*sum([α[a][i]*α[b][j]*α[a][k] + α[b][i]*α[a][j]*α[b][k]
                         for i = 1:N-1 for j = 1:N for k = i+1:N])
               ) +
               16/N^4 * (
                    sum([α[a][i]^2*α[b][j]^2 for i = 1:N for j = 1:N]) +
                    2*sum([α[a][i]^2*α[b][j]*α[b][k] + α[b][i]^2*α[a][j]*α[a][k]
                           for i = 1:N for j = 1:N-1 for k = j+1:N]) +
                    4*sum([α[a][i]*α[b][j]*α[a][k]*α[b][l]
                           for i = 1:N-1 for j = 1:N-1
                           for k = i+1:N for l = j+1:N])
               )
        )
    end
end

#Assert that equation A18 is identical to just ra*rb
#(except for rounding error due to floating-point arithmetic).
@assert maximum(abs.(rr_all - rarb)) < 1e-12

#αsums store the final sum in equation A18
αsums = []
for a = 1:Nₐ
    for b = 1:Nₐ
        push!(αsums, sum([α[a][i]*α[b][j]*α[a][k]*α[b][l]
              for i = 1:N-1 for j = 1:N-1 for k = i+1:N for l = j+1:N]))
    end
end

#Left-hand side of equation A22
Eu_αsums_LHS_eq_A22 = Dict(
    [(u, sum([αsums[j] for j = 1:Nₐ^2 if u_all[j] == u]) / nᵤ[u]) for u=0:N2])

#Equation A19
VarMSE𝒦ₕt_eqA19 = 64 / N^4 * sum(
    [(vᵤ𝒦ₕt[u] - nᵤ[u] / Nₐ^2) * Eu_αsums_LHS_eq_A22[u] for u = 0:N2])

#Generate the sums for equation A22
αsums_iijj = []
αsums_iijk = []
αsums_ijjk = []
αsums_ijkk = []
αsums_ijkl = []
αsums_ikjl = []
αsums_iljk = []
for a = 1:Nₐ
    for b = 1:Nₐ
        push!(αsums_iijj, sum([α[a][i]*α[b][i]*α[a][j]*α[b][j]
                               for i = 1:N-1 for j = i+1:N]))
        push!(αsums_iijk, sum([α[a][i]*α[b][i]*α[a][j]*α[b][k]
                               for i = 1:N-2 for j = i+1:N-1 for k=j+1:N]))
        push!(αsums_ijjk, sum([α[a][i]*α[b][j]*α[a][j]*α[b][k]
                               for i = 1:N-2 for j = i+1:N-1 for k=j+1:N]))
        push!(αsums_ijkk, sum([α[a][i]*α[b][j]*α[a][k]*α[b][k]
                               for i = 1:N-2 for j = i+1:N-1 for k=j+1:N]))
        push!(αsums_ijkl, sum([α[a][i]*α[b][k]*α[a][j]*α[b][l]
              for i = 1:N-3 for j = i+1:N-2 for k=j+1:N-1 for l=k+1:N]))
        push!(αsums_ikjl, sum([α[a][i]*α[b][j]*α[a][k]*α[b][l]
              for i = 1:N-3 for j = i+1:N-2 for k=j+1:N-1 for l=k+1:N]))
        push!(αsums_iljk, sum([α[a][i]*α[b][j]*α[a][l]*α[b][k]
              for i = 1:N-3 for j = i+1:N-2 for k=j+1:N-1 for l=k+1:N]))
    end
end

#Right-hand side of equation A22
Eu_αsums_RHS_eq_A22 = Dict(
    [(u, sum([αsums_iijj[j] + 2*αsums_iijk[j] + 2*αsums_ijjk[j] +
     2*αsums_ijkk[j] + 2*αsums_ijkl[j] + 2*αsums_ikjl[j] + 2*αsums_iljk[j]
     for j = 1:Nₐ^2 if u_all[j] == u]) / nᵤ[u]) for u=0:N2])

#Assert that equation that the left-hand side in equation A22 is
#the same that the right-hand side (except for rounding error
#due to floating-point arithmetic).
for u = 0:N2
    @assert abs(Eu_αsums_LHS_eq_A22[u] - Eu_αsums_RHS_eq_A22[u]) < 1e-12
end

#To create the sums below, it is useful to store the potential outcomes
#in a Nx2 matrix called Y
Y = hcat(Y0, Y1)

#Generate the sums in equations A30, A31 and A32
S₁ = sum([Y[i, x₁+1]*Y[j, x₂+1]*Y[i, x₃+1]*Y[j, x₄+1]
          for x₁=0:1 for x₂=0:1 for x₃=0:1 for x₄=0:1 for i=1:N-1 for j=i+1:N])
S₂ = sum([Y[i, x₁+1]*Y[j, x₂+1]*Y[i, x₃+1]*Y[k, x₄+1] +
          Y[i, x₁+1]*Y[j, x₂+1]*Y[j, x₃+1]*Y[k, x₄+1] +
          Y[i, x₁+1]*Y[k, x₂+1]*Y[j, x₃+1]*Y[k, x₄+1]
          for x₁=0:1 for x₂=0:1 for x₃=0:1 for x₄=0:1
          for i=1:N-2 for j=i+1:N-1 for k=j+1:N])
S₃ = sum([Y[i, x₁+1]*Y[j, x₂+1]*Y[k, x₃+1]*Y[l, x₄+1] +
          Y[i, x₁+1]*Y[k, x₂+1]*Y[j, x₃+1]*Y[l, x₄+1] +
          Y[i, x₁+1]*Y[l, x₂+1]*Y[j, x₃+1]*Y[k, x₄+1]
          for x₁=0:1 for x₂=0:1 for x₃=0:1 for x₄=0:1
          for i=1:N-3 for j=i+1:N-2 for k=j+1:N-1 for l=k+1:N])

ᴺP₂ = N*(N-1)
ᴺP₃ = N*(N-1)*(N-2)
ᴺP₄ = N*(N-1)*(N-2)*(N-3)
ψ = S₁ / ᴺP₂ - 2 * S₂ / ᴺP₃ + 4 * S₃ / ᴺP₄

#Equation A34
VarMSE𝒦ₕt_eqA34 = 64 / N^4 * ψ * sum([(vᵤ𝒦ₕt[u] - nᵤ[u] / Nₐ^2) * (u-N/4)^2
                                      for u = 0:N2])

#Equation A35
VarMSE𝒦ₕt_eqA35 = 64 / N^4 * ψ * ((Nₐ-H) / (H*Nₐ) * N^2/8 +
                              sum([(vᵤ𝒦ₕt[u] - nᵤ[u] / Nₐ^2) * (u-N/4)^2
                                   for u = 1:N2-1]))

#Equation 14
ϕ𝒦ₕt = (4/N)^2 *H/(H-2) *  sum([vᵤ𝒦ₕt[u] * (u-N/4)^2 for u=1:N2-1])
#Equation A36
ϕ𝒦 = (4/N)^2 * Nₐ/(Nₐ-2) * sum([nᵤ[u] / Nₐ^2 * (u-N/4)^2 for u=1:N2-1])

#Theorem 4 (equation 13)
VarMSE𝒦ₕt_eq13 = 4/N^2*ψ*(2*(Nₐ-H)/(H*Nₐ)+ (H-2)/H * ϕ𝒦ₕt - (Nₐ-2)/Nₐ *ϕ𝒦)

println("********************")
println("Print the variances of the MSE for the set 𝒦ₕ from different",
        "equations.\nThese should all be identical")
println(VarMSE𝒦ₕ)
println(VarMSE𝒦ₕ_eqA8)
println(VarMSE𝒦ₕ_eqA9)
println(VarMSE𝒦ₕ_eqA11)
println("********************")
println("Print the variances of the MSE for the set 𝒦ₕt from different",
        "equations\nThese should all be identical")
println(VarMSE𝒦ₕt)
println(VarMSE𝒦ₕt_eqA14)
println(VarMSE𝒦ₕt_eqA15)
println(VarMSE𝒦ₕt_eqA19)
println(VarMSE𝒦ₕt_eqA34)
println(VarMSE𝒦ₕt_eqA35)
println(VarMSE𝒦ₕt_eq13)
