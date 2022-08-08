#This code has been tested with Julia version 1.7.2
using DataFrames #Tested with version 1.3.4
using CSV #Tested with version 0.10.4
using Statistics
using StatsPlots #Tested with version v0.14.34
using Distributions #Tested with version v0.23.12
using GLM #Tested with version 1.8.0
using OrderedCollections #Tested with version 1.4.1

function gen_exp_phi(n)
    #Calculate φ under complete randomization
    na = binomial(BigInt(n), BigInt(n/2))
    n2 = Int(n/2)
    return Float64(sum([
        binomial(n2,u)*binomial(n2,n2-u)/na * (4/n)^2 * (u - n/4)^2
        for u = 1:n2-1]))
end

function get_varmse(psi, phi_d, n, h)
    #Calculate the varianec of the MSE under the veil of ignorance
    na = binomial(BigInt(n),BigInt(n/2))
    varvar = 4 / n^2 * psi *(2 * (na - h) / h / na + (h-2)/h * phi_d -
          (na - 2)/ na * gen_exp_phi(n))
    return varvar
end

function get_violin_plot(distr, rsq, f, xt, jl, dsgdict, lg, xl, yl, yt, el)
    #Make the violin plots (figures 2 and 3 in the paper)
    dft = dfm[dfm[!, :distr] .== distr, :]
    xtickvals = [x * f for x = xt]
    StatsPlots.violin()
    x = Float64[]
    y = Float64[]
    for (dsg, dsgt) = dsgdict
        push!(x, mean(dft[!, string("as_cor", dsgt)]))
        push!(y, mean(dft[!, string("std_mse_", dsgt, "__", RSQDICT[rsq])]))
        StatsPlots.violin!([mean(dft[!, string("as_cor", dsgt)]) * f],
                           dft[!, string("std_mse_", dsgt, "__", RSQDICT[rsq])],
                           label=string(dsg))
        if jl != false
            jitter = rand(Uniform(-jl, jl), size(dft)[1])
            plot!([mean(dft[!, string("as_cor", dsgt)]) * f] .+ jitter,
                  dft[!, string("std_mse_", dsgt, "__", RSQDICT[rsq])],
                  label = "", seriestype = :scatter, markersize=2,
                  markeralpha=.5, markerstrokewidth=0, markercolor="gray")
        end
    end
    if lg == true
        plot!(legend=:outertopright, xlabel = "", ylabel="", legendtitle="pa",)
    else
        plot!(legend=false, xlabel = "", ylabel="")
    end
    if xt == false
        plot!(xticks = false)
    else
        plot!(xticks = (xtickvals, [x / f for x = xtickvals]))
    end
    if yt == false
        plot!(yticks = false)
    end
    if el == true
        plot!(x .* f, y, color="gray", label=false)
        scatter!(x .* f, y, color="black", label=false)
    end
    plot!(xlims = (xl[1] * f, xl[2] * f))
    plot!([xl[1] * f, xl[2] * f], [1-rsq, 1-rsq], linestyle=:dash, color="gray",
          label=false)
    plot!(ylims = (yl[1], yl[2]))
    plot!(ylims = (yl[1], yl[2]))
    plot!(title=string("R-squared = ", rsq), titlefontsize=10, grid=false)
end

function eval_var_plot(rsq, distr)
    #Comparison between the theoretical and empirical variances of the MSE for
    #the simulation study in Section 3.5 in the paper. Generates figure A2
    #in the supplementary material.
    p = Plots.plot()
    A = []
    B = []
    dft = dfm[dfm[!, :distr] .== distr, :]
    psi_m = 8
    for (dsg, dsgt) = DSGDICT
        m = mean(dft[!, string("as_cor", dsgt)])
        q1 = var(dft[!, string("std_mse_", dsgt, "__", RSQDICT[rsq])])
        q2 = get_varmse(mean(
            dft[!, string("psi__", RSQDICT[rsq])]), m, N, H) / exp_CR_mse^2
        push!(A, q1)
        push!(B, q2)
    end
    B2 = deepcopy(B)
    push!(B2, 0)
    p = Plots.plot!(B, A, seriestype = :scatter, markerstrokewidth=1.8)
    p = Plots.plot!(B2, B2)

    plot!(xlims = (0, maximum(B) * 1.05),
          ylims = (0, maximum(B) * 1.05))
    a = plot!(title = "", legend=false, grid=false, size = (350, 300))
    return p
end

H = 10000
rsq = 0.0
N = 50
K = 5
exp_CR_mse = 4 / N

PALIST = [1, .5, .25, .1, .05, .025, .01, .005, .0025, .001]
RSQLIST = 0:.05:.15
DSGLIST = vcat(PALIST, "ps-alg")

DSGDICT = OrderedDict(
    [(des, replace(string(des), "." => "_")) for des = DSGLIST])
DSGDICT["ps-alg"] = "alg"
RSQDICT = OrderedDict(
    [(rsq, replace(string(rsq), "." => "_")) for rsq = RSQLIST])

infile_normal = string("sim_var_mse_fixed_x_normal.csv")
infile_lognormal = string("sim_var_mse_fixed_x_lognormal.csv")

dfnormal = CSV.read(infile_normal, DataFrame)
dfnormal[!, :distr] .= "Normal"
dflognormal = CSV.read(infile_lognormal, DataFrame)
dflognormal[!, :distr] .= "Lognormal"
dfm = vcat(dfnormal, dflognormal)
dfm = dfm[dfm[!, :H] .== H, :]

mse_list_r2_list = vcat([[string("mse_", dsgt, "__", rsqtext)
                         for (dsg, dsgt) = DSGDICT]
                         for (i, rsqtext) = RSQDICT]...)
as_cor_list = [string("as_cor", dsgt) for (dsg, dsgt) = DSGDICT]

for i = 1:length(mse_list_r2_list)
    dfm[!, string("std_", mse_list_r2_list[i])] .=
    dfm[!, mse_list_r2_list[i]] ./ exp_CR_mse
end

s1 = get_violin_plot("Normal", 0, 6000, .02:.001:.023, false,
    DSGDICT, false, (0.02, 0.023), (0.3, 1.2), true, true)
s2 = get_violin_plot("Normal", 0.05, 6000, .02:.001:.023, false,
    DSGDICT, true, (0.02, 0.023),(0.3, 1.2), true, true)
s3 = get_violin_plot("Normal", 0.1, 6000, .02:.001:.023, false,
    DSGDICT, false, (0.02, 0.023), (0.3, 1.2), true, true)
s4 = get_violin_plot("Normal", 0.15, 6000, .02:.001:.023, false,
    DSGDICT, false, (0.02, 0.023), (0.3, 1.2), true, true)

a = plot(s1, s2, s3, s4, layout = @layout([A{.44w} B{.56w} ; C{.44w} D{.44w} _]),
     size = (700, 600))

savefig(string("figure2.pdf"))

xl = (0.02, 0.11)
xlist = .025:.025:.1
yl = (0, 5)
f = 500

s1 = get_violin_plot("Lognormal", 0, f, xlist, false,
    DSGDICT, false, xl, yl, true, true)
s2 = get_violin_plot("Lognormal", 0.05, f, xlist, false,
    DSGDICT, true, xl,yl, true, true)
s3 = get_violin_plot("Lognormal", 0.1, f, xlist, false,
    DSGDICT, false, xl, yl, true, true)
s4 = get_violin_plot("Lognormal", 0.15, f, xlist, false,
    DSGDICT, false, xl, yl, true, true)

plot(s1, s2, s3, s4, layout = @layout([A{.44w} B{.56w} ; C{.44w} D{.44w} _]),
          size = (700, 600))

savefig(string("figure3.pdf"))

eval_var_plot(0, "Normal")
savefig(string("figureA2a.pdf"))

eval_var_plot(0, "Lognormal")
savefig(string("figureA2b.pdf"))
