#This code has been tested with Julia version 1.7.2
using Random
using DataFrames #Tested with version 1.3.4
using CSV #Tested with version 0.10.4
using Dates
using GLM #Tested with version 1.8.0
using Statistics
using LinearAlgebra

function get_as_cor(wset)
    #Calculate the assignment correlation
    n = length(wset[1])
    sqop = 0
    h2 = length(wset)
    n4 = n / 4
    for i = 1:h2 - 1
        for j = i+1:h2
            biop = sum(abs.(wset[i] .- wset[j]))
            sqop += (n4 - biop/2)^2
        end
    end
    avg_sqop = sqop / binomial(h2, 2)
    return avg_sqop * 16 / n^2
end

function print_time(rep_nr, tot_reps)
    #Output information on how many simulations that have been performed
    #and when the simulation is expected to finish
    elapsed_time = time_ns() - START_TIME
    perc_compl = (rep_nr / tot_reps) * 100
    exp_tot_time = elapsed_time / (perc_compl / 100)
    eft = Dates.now() + Dates.Nanosecond(floor(exp_tot_time - elapsed_time))
    perc_compl = round(perc_compl, sigdigits=3)
    println("", perc_compl, "% completed, expect to finish ", eft)
end

function calc_mahal(x_t, x_mean, vi, n, ones)
    #Function to calculate the Mahalanobis distance
    x_t_mean = transpose(x_t) * ones
    diff_x = 2 * (x_t_mean - x_mean)
    m = dot(diff_x, vi * diff_x)
    return n / 4 * m
end

function base_alg(x0, x1, x0_ind, x1_ind, vi, n, x_mean, ones)
    #Pair-switching algorithm of Krieger et al (2019)
    m_b = calc_mahal(x1, x_mean, vi, n, ones)
    n2 = Int(n / 2)
    switch_made = 1
    nr_calc = 0
    local y0_ind
    local y1_ind
    local y0
    local y1
    while switch_made == 1
        switch_made = 0
        for i = 1:n2
            for j = 1:n2
                sl1 = x1[j, :]
                sl0 = x0[i, :]
                x0[i, :] = sl1
                x1[j, :] = sl0
                m = calc_mahal(x1, x_mean, vi, n, ones)
                nr_calc += 1
                if m < m_b
                    switch_made = 1
                    m_b = m
                    y0_ind = copy(x0_ind)
                    y1_ind = copy(x1_ind)
                    tl1 = y1_ind[j]
                    tl0 = y0_ind[i]
                    y0_ind[i] = tl1
                    y1_ind[j] = tl0
                    y0 = copy(x0)
                    y1 = copy(x1)
                end
                x0[i, :] = sl0
                x1[j, :] = sl1
            end
        end
        if switch_made == 1
            x0_ind, x1_ind = y0_ind, y1_ind
            x0, x1 = y0, y1
        end
    end
    return x0, x1, nr_calc, m_b, x0_ind, x1_ind
end

function rerand_alg(x, h)
    #Get h assignment vectors with the pair-switching algorithm of
    #Krieger et al. (2019).
    n, k = size(x)
    h2 = Int(h / 2)
    n2 = Int(n / 2)
    vi = inv(cov(x))
    ones_full = (zeros(n) .+ 1) ./ n
    ones_half = (zeros(n2) .+ 1) ./ n2
    x_mean = transpose(x) * ones_full
    t1 = Array{Int64, 1}(undef, n2)
    fill!(t1, 0)
    WLIST = Array{Array{Int64, 1}}(undef, h2)
    fill!(WLIST, t1)
    M = Array{Float64}(undef, h2)
    nr_found = 0
    nr_tried = 0
    max_tries = h2 * 100
    while (nr_found < h2)
        if nr_tried == max_tries
            return WLIST[1:nr_found], nr_tried, maximum(M[1:nr_found])
        end
        rand_samp = shuffle(2:n)
        x1_ind = append!([1], rand_samp[1:n2-1])
        x0_ind = rand_samp[n2:n-1]
        x1_start = x[x1_ind, :]
        x0_start = x[x0_ind, :]
        q = base_alg(x0_start, x1_start, x0_ind, x1_ind,
                      vi, n, x_mean, ones_half)
        if 1 in q[5]
            t_ind = copy(q[5])
        elseif 1 in q[6]
            t_ind = copy(q[6])
        else
            throw(ValueError)
        end
        nr_tried += 1
        W = Int.(zeros(n))
        W[t_ind] .+= 1
        if (W in WLIST) == false
            nr_found += 1
            M[nr_found] = q[4]
            WLIST[nr_found] = W
        end
    end
    return WLIST, nr_tried, maximum(M), M
end

function rerand(x, h, palist)
    #Get h assignment vectors for the rerandomization designs, given a
    #list of values for pₐ
    n, k = size(x)
    vi = inv(cov(x))
    h2 = Int(h / 2)
    n2 = Int(n / 2)
    ones_full = (zeros(n) .+ 1) ./ n
    ones_half = (zeros(n2) .+ 1) ./ n2
    x_mean = transpose(x) * ones_full
    nr_rerand_list = [Int.(h2 / pa) for pa = palist]
    nr_rerand = maximum(nr_rerand_list)
    M = Array{Float64}(undef, nr_rerand)
    T_IND = Array{Array}(undef, nr_rerand)
    for i = 1:nr_rerand
        rand_samp = shuffle(2:n)
        t_ind = append!([1], rand_samp[1:n2-1])
        xt = x[t_ind, :]
        M[i] = calc_mahal(xt, x_mean, vi, n, ones_half)
        T_IND[i] = t_ind
    end
    WLIST_ALL = Array[]
    M_MAX_ALL = Float64[]
    MLIST = Array[]
    for i = nr_rerand_list
        c_r = sortperm(M[1:i])[1:h2]
        push!(MLIST, M[c_r])
        T_IND_U = T_IND[c_r]
        WLIST = Array[]
        for j = T_IND_U
            t = Int.(zeros(n))
            t[j] .+= 1
            push!(WLIST, t)
        end
        push!(WLIST_ALL, WLIST)
        push!(M_MAX_ALL, maximum(M[c_r]))
    end
    return unique(WLIST_ALL), nr_rerand_list, M_MAX_ALL, MLIST
end

function mse_set(y, wset)
    #Calculate MSE for a set of assignment vectors
    return mean([(mean(y[i .== 1]) - mean(y[i .== 0]))^2 for i = wset])
end

function sim_study(x, rsqlist, h, palist)
    #Generate the set of assignment vectors for each design
    wset_rerand_set, nr_tried_rerand_set, m_max_rerand_set, m_list_rerand =
        rerand(x, h, palist)
    wset_alg, nr_tried_alg, m_max_alg, m_list_alg = rerand_alg(x, h)
    wsetset = vcat(wset_rerand_set, [wset_alg])
    nrtried = vcat(nr_tried_rerand_set, [nr_tried_alg])
    mmax = vcat(m_max_rerand_set, [m_max_alg])
    mlist = vcat(m_list_rerand, [m_list_alg])
    hset = [length(wset)*2 for wset = wsetset]
    ascorlist = [get_as_cor(wset) for wset = wsetset]
    return vcat(nrtried, mmax, hset, ascorlist), wsetset
end

function gen_namelist(palist, rsqlist)
    #Create variable names for the dataframe where results are stored
    endlist = vcat([replace(string(i), "." => "_") for i = palist], ["alg"])
    rsqtextlist = [replace(string(rsq), "." => "_") for rsq = rsqlist]
    stublist = ["nr_tried_", "m_max_", "act_h", "as_cor"]
    out = [j * i for j = stublist for i = endlist]
    mse_list = [string("mse_", i, "__", r) for i = endlist for r = rsqtextlist]
    return vcat(out, mse_list)
end

function get_data(file, row_nr; matrix=true)
    #Read the data from the file generated by sim_rerand_alg.jl.
    #Choose given row number. To replicate the findings in the paper,
    #row_nr 1 should be chosen for both normal and lognormal covariates.
    f = readlines(file)[row_nr]
    f = replace(f, "[" => "")
    f = replace(f, "]" => "")
    if matrix == true
        f = split(f, "; ")
        s = [parse.(Float64, split(i, " ")) for i = f]
        return Matrix(transpose(hcat(s...)))
    else
        s = parse.(Float64, split(f, ", "))
        return reshape(s, length(s), 1)
    end
end

function gen_data(x, e, rsqlist, vare)
    ylist = []
    n, k = size(x)
    beta = reshape(ones(k), k, 1)
    for rsq = rsqlist
        if rsq == 0.0
            y = e
        else
            c = sqrt(var(x*beta) * (1-rsq) / vare / rsq)
            y = x * beta + c * e
            y = reshape(y, n)
        end
        y = y ./ sqrt(var(y))
        push!(ylist, y)
    end
    return ylist
end

function get_psi_homo(Y, N)
    S₃ = sum([Y[i]^2*Y[j]^2 for i=1:N-1 for j=i+1:N])
    S₄ = sum([Y[i]^2*Y[j]*Y[k] + Y[i]*Y[j]^2*Y[k] + Y[i]*Y[j]*Y[k]^2
              for i=1:N-2 for j=i+1:N-1 for k=j+1:N])
    S₅ = sum([Y[i]*Y[j]*Y[k]*Y[l]
             for i=1:N-3 for j=i+1:N-2 for k=j+1:N-1 for l=k+1:N])
    S̄₃ = S₃ / binomial(N, 2)
    S̄₄ = S₄ / (3 * binomial(N, 3))
    S̄₅ = S₅ / binomial(N, 4)
    ψ = 8 * (S̄₃ - 2 * S̄₄ + S̄₅)
    return ψ
end

function get_data_row(distr, n)
    #Return the row where the data used in the simulations are stored.
    #For the normal distribution, the first set of covariates generated
    #with the seed 12345 in the script sim_rerand_alg.jl is used,
    #whereas for the lognormal distribution, the first set of covariates
    #generated with the seed 1234 in the script sim_rerand_alg.jl is used.
    if distr == "normal"
        gen_seed = 12345
    elseif distr == "lognormal"
        gen_seed = 1234
    end
    inpfile = string("rerand_alg_", distr, "_", string(n), ".csv")
    df = CSV.read(inpfile, DataFrame)
    df[!, :rows] = 1:10000
    row_nr = df[df[!, :seed] .== gen_seed, :rows][1]
    return row_nr
end

cd("D:\\Dropbox\\NN design\\tex\\PAPER\\Properties-of-restricted-randomization")

PALIST = [1, .5, .25, .1, .05, .025, .01, .005, .0025, .001]
RSQLIST = 0:.05:.15
namelist = gen_namelist(PALIST, RSQLIST)
rsqtextlist = [replace(string(rsq), "." => "_") for rsq = RSQLIST]

H = 10000
REPS = 10000
N = 50
σe = 1
seed = 12345
START_TIME = 0
for distr = ["normal", "lognormal"]
    fileX = string("rerand_alg_", distr, "_X_", string(N), ".csv")
    outfile = string("sim_var_mse_fixed_x_", distr, ".csv")
    data_row = get_data_row(distr, N)
    Random.seed!(seed)
    X = get_data(fileX, data_row, matrix=true)
    q = sim_study(X, RSQLIST, H, PALIST)
    global START_TIME = time_ns()
    for i = 1:REPS
        e = rand(Normal(0, σe), N)
        YLIST = gen_data(X, e, RSQLIST, σe^2)
        mse_list = [mse_set(y, wset) for wset = q[2] for y = YLIST]
        qn = vcat(q[1], mse_list)
        df = DataFrame(reshape(qn, 1, length(qn)), :auto)
        rename!(df, Symbol.(namelist))
        for j = 1:length(RSQLIST)
            ols = GLM.fit(LinearModel, [ones(N) X], YLIST[j])
            df[!, string("rsq__", rsqtextlist[j])] .= r2(ols)
            df[!, string("psi__", rsqtextlist[j])] .= get_psi_homo(YLIST[j], N)
        end
        df[!, :seed] .= seed
        df[!, :data_row] .= data_row
        df[!, :H] .= H
        if isfile(outfile) == false
            CSV.write(outfile, df)
        else
            CSV.write(outfile, df, append=true)
        end
        print_time(i, REPS)
    end
end
