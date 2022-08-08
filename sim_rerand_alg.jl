#This code has been tested with Julia version 1.5.1. Note that later versions
#of Julia will have a different random number generator
using Distributions #Tested with version 0.23.8
using Random
using LinearAlgebra
using DataFrames #Tested with version 1.2.2
using CSV #Tested with version 0.9.9
using Dates

function restart(file, seed, a, n)
    #The simulations take a long time. If for some reason, the simulation
    #is aborted before finishing, this function allows a restart for a
    #given position with the given seed
    df = CSV.read(file, DataFrame)
    dfs = df[df[!, :seed] .== seed, [:nr_rerand4, :nr_alg]]
    dfs[!, :tot] = dfs[!, :nr_rerand4] .+ dfs[!, :nr_alg]

    for i = dfs[!, :tot]
        X = transpose(rand(a, n))
        for j = 1:i
            shuffle(2:n)
        end
    end
    return size(dfs)[1]
end

function get_as_cor(wset, n)
    #Calculate the assignment correlation
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
    println(perc_compl, "% completed, expect to finish ", eft)
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
    return WLIST, nr_tried, maximum(M)
end

function rerand(x, h, pa_list)
    #Get h assignment vectors for the rerandomization designs, given a
    #list of values for pₐ
    n, k = size(x)
    vi = inv(cov(x))
    h2 = Int(h / 2)
    n2 = Int(n / 2)
    ones_full = (zeros(n) .+ 1) ./ n
    ones_half = (zeros(n2) .+ 1) ./ n2
    x_mean = transpose(x) * ones_full
    nr_rerand_list = [Int.(h2 / pa) for pa = pa_list]
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
    for i = nr_rerand_list
        c_r = sortperm(M[1:i])[1:h2]
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
    return unique(WLIST_ALL), nr_rerand_list, M_MAX_ALL
end

function sim_study()
    #Perform the simulation study
    reps_remaining = TOT_REPS - nr_done
    fl = Int.(floor(reps_remaining / inner_rep_target))
    rem = reps_remaining%inner_rep_target
    INNER_REPS = Int.(zeros(fl) .+ inner_rep_target)
    INNER_REPS[1:rem] .+= 1
    c = 0
    for i = INNER_REPS
        OUT = Array{Any}(undef, i, length(namelist))
        q = ""
        for k = 1:i
            X = transpose(rand(A, N))
            wset_rerand_set, nr_tried_rerand_set, m_max_rerand_set =
                rerand(X, H, PA_LIST)
            out = Any[]
            for j = 1:length(PA_LIST)
                append!(out, get_as_cor(wset_rerand_set[j], N))
                append!(out, nr_tried_rerand_set[j])
                append!(out, m_max_rerand_set[j])
                append!(out, length(wset_rerand_set[j])*2)
            end
            wset_alg, nr_tried_alg, m_max_alg = rerand_alg(X, H)
            as_cor_alg = get_as_cor(wset_alg, N)
            h_alg = length(wset_alg)*2
            append!(out, [as_cor_alg, nr_tried_alg, m_max_alg, h_alg, seed])
            OUT[k, :] = out
            q = q * string(X) * "\n"
        end
        df = DataFrame(OUT, :auto)
        rename!(df, Symbol.(namelist))
        if isfile(FILENAME) == false
            CSV.write(FILENAME, df)
        else
            CSV.write(FILENAME, df, append=true)
        end
        io2 = open(FILENAME2, "a")
        write(io2, q)
        close(io2)
        c += i
        print_time(c, sum(INNER_REPS))
    end
end

#Because the script takes a long time to run, to efficiently run in parallel
#while still make results reproducible, 5 different scripts was run with
#5 different seed values with 2000 reps each (TOT_REPS constant below).
#The seed values chosen were 1, 12, 123, 1234 and 12345.
seed = 1234
N = 50
H = 10000
K = 5
PA_LIST = [1, 0.1, 0.01, 0.001]
TOT_REPS = 2000
inner_rep_target = 1

MV = zeros(K)
CM = Matrix{Float64}(I, K, K)

#Set distribution (normal or lognormal)
distr = "lognormal"
if distr == "normal"
    A = MvNormal(MV, CM)
elseif distr == "lognormal"
    A = MvLogNormal(MV, CM)
end

namelist = []
for i = 1:length(PA_LIST)
    append!(namelist, ["as_cor_rerand" * string(i)])
    append!(namelist, ["nr_rerand" * string(i)])
    append!(namelist, ["m_max_rerand" * string(i)])
    append!(namelist, ["act_h" * string(i)])
end

append!(namelist, ["as_cor_alg", "nr_alg", "m_max_alg", "act_h_alg", "seed"])

FILENAME = string("rerand_alg_", distr, "_", string(N), ".csv")
FILENAME2 = string("rerand_alg_", distr, "_X_", string(N), ".csv")

START_TIME = time_ns()
Random.seed!(seed)

#Restart if necessary (see restart function)
if isfile(FILENAME) == true
    nr_done = restart(FILENAME, seed, A, N)
else
    nr_done = 0
end

sim_study()
