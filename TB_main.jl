
#module TB
using Parameters, Distributions, StatsBase, StaticArrays, Random, Match, DataFrames, DelimitedFiles, CSV

#health status of the individuals
@enum HEALTH SUSC = 0 LAT = 1 LATT = 2 ACT = 3  ACTT = 4 DS = 5 UNDEF = 6
#LATT = Latent under treatment , ACTT = Active under treatment , DS = Dormant State

#=-----------------------------
Construct each population member
------------------------------=#

Base.@kwdef mutable struct Human

    id::Int64 = 0 #We have 39353 person
    AgeGroup::Int16 = 0
    Age::Int16 = 0 
    AGC::Int16 = 0 #Age Group for Contact Matrix
    AGR::Int16 = 0 #Age Group for last Report

    Health:: HEALTH = SUSC
    swap::HEALTH = UNDEF #What will be their next state
    sickfrom::Int64 = 0
    tis::Int64 = 0 #Time in the current state
    exp::Int64 = 0 #Max time in the state
    α:: Float64 = 0
    nextweek_meets:: Int16 = 0
    Hid::Int64 = 0

    vaccinated::Bool = false
    treated::Bool = false
    act_hist::Bool = false
    will_act::Bool = false
    iso :: Bool = false
    rep_contact :: Array{Int64} = [] #This is the array for repeated contacts
    same_hh :: Array{Int64} = [] #This is the array for same household members of each individual
    HIV::Bool = false
    Diab::Bool = false
    renal::Bool = false

end


#=------------------------------
----Construct each Household----
------------------------------=#

Base.@kwdef mutable struct Household
    Hid:: Int64 = 0 # We have 10356 household
    Size:: Int64 = 0
    AvailableSize:: Int64 = 0
    Infected::Bool = false
    # Adult::Bool #Each household should have at least one adult
end

@with_kw mutable struct ModelParameters @deftype Float64
    # general parameters
    sim_time::Int64 = 7800  ## 110 years: 70 warmup, 10 calibration, 30 Model run

    vaccine_onoff::Bool = false
    treatment_onoff::Bool = false
    calibration:: Bool = false
    test_sens :: Float64 = 0
    compliance :: Float64 = 0
    efficacy :: Float64 = 0

    β = 0.31 ## We get this by calibration
    θ = 0.32    ## Relapse rate (this can change by the history of the person)
    α_red = 0 #the reducion in activation probability of the person based on treatment scenario
    trace_cov = 0
    vaccin_cov = 0   ##Vaccine coverage (will change in different scenarios)
    lat_treat_time = 0  ## Different scenarios? #3 6 and 9 months
    act_treat_time = 0 ## Different scenarios? #4 or 7 months (Regimen 1 or 2)
    initial_infected::Int64 = 5
    diagnosis :: Int64 = 4

end


#=------------------------------x
--------global variables--------
------------------------------=#

const total_population = 39353 #According to gov of Nunavut
const total_household  = 10356 #According to statista
const Age_dist = Categorical(@SVector[0.115732369,0.114619558,0.09500626,0.087355682,0.084712756,0.088607595,0.078314091,0.06635137,
0.05647517,0.05744888,0.051189317,0.039087495,0.027124774,0.018222284,0.009597997,0.005703158,0.00292113,0.001112811,0.000278203,0.000139101])
const Age_braks= @SVector [0:4,5:9,10:14,15:19,20:24,25:29,30:34,35:39,40:44,45:49,50:54,55:59,60:64,65:69,70:74,75:79,80:84,85:89,90:94,95:99]
const latent_dist = Categorical(@SVector[0.17188727,0.15832081,0.08204884,0.067979255,0.07089221,0.061512495,0.080964945,0.07997762,0.10417054,0.056408,0.065838015])
const latent_braks = @SVector[1:15, 16:30, 31:45, 56:60, 61:80, 81:100, 101:150, 151:200, 201:300, 301:400, 401:5000]
const human = Array{Human}(undef,total_population)
const household = Array{Household}(undef,total_household)
const p = ModelParameters()
export ModelParameters, HEALTH, Human, human

#-----------Init
function init_human()
    @inbounds for i in 1:total_population
        human[i]=Human()
        x = human[i]
        x.id = i
        x.exp = 99999
        #get_nextday_counts

    end
end

function init_household()
    @inbounds for i in 1:total_household
        household[i] = Household()
        household[i].Hid = i
        # household[i].Adult = 0
    end
    Smin = [1,    1952, 3945, 5522, 7157]
    Smax = [1951, 3944, 5521, 7156, 10356]

    sum_population = 0

    for H in 1:5
        sum_population += H*(Smax[H]-Smin[H]+1)


        for j in Smin[H]:Smax[H]
            household[j].Size = H
            household[j].AvailableSize=H
        end

        #This is for the last group size , in previous loop we have
        #already assigned 5 people to each household.

        if H == 5
            remain_pop = length(human)-sum_population
            while remain_pop !== 0
                if remain_pop >= 15
                    remain_size = rand(0:15)
                else
                    remain_size = rand(0:remain_pop)
                end
                # while length(filter(x -> x.Size <= 20, household)) !== 0
                remain_id = rand(Smin[5]:Smax[5])
                if household[remain_id].Size >= 6
                    household[remain_id].Size +=0
                    household[remain_id].AvailableSize += 0
                    remain_pop += 0
                else
                    household[remain_id].Size += remain_size
                    household[remain_id].AvailableSize += remain_size
                    remain_pop -= remain_size
                end
            end
        end
    end
end


#--------Apply Age groups to each Agent-------#

function apply_age_g(x)
    agegp=rand(Age_dist)::Int64
    x.AgeGroup=agegp
    x.Age = rand(Age_braks[x.AgeGroup])
    x.AGC = get_AGC(x)
    x.AGR = get_AGR(x)
    x.α   = get_alpha(x)
    get_comorb(x)
    get_alpha_comorb(x)
end
export apply_age_g

function get_AGC(x)
@match (x.Age) begin
    0:4     => 1
    5:19    => 2
    20:49   => 3
    50:64   => 4
    64:99   => 5
    _       => println("No Age Group Specified")
end
end

function get_AGR(x)
@match (x.Age) begin
    0:15     => 1
    16:99    => 2
    _       => println("No Age Group Specified")
end
end

function get_comorb(x::Human)
x.Diab = rand()< 0.04 ? true : false
x.HIV = rand()< 0.0025 ? true : false
x.renal = rand() < 0.01511 ? true : false
end

function get_alpha(x)
    @match (x.Age) begin
        0:1   => rand(Uniform(0.3,0.4))
        1:2   => rand(Uniform(0.1,0.2))
        2:5   => 0.05
        5:10  => 0.02
        10:14 => rand(Uniform(0.1,0.2))
        15:25 => 0.13
        26:35 => 0.12
        36:45 => 0.07
        46:55 => 0.06
        56:65 => 0.03
        66:99 => rand(Uniform(0.05,0.1))
    end
end


function get_alpha_comorb(x)
if x.HIV == true
    x.α = (9.9)*x.α
    else
    x.α = x.α
end

if x.Diab == true
    x.α = (1.7)*x.α
    else
    x.α = x.α
end

if x.renal == true
    x.α = (2.4)*x.α
    else
    x.α = x.α
end
if x.α > 1
x.α = 1
end
end

## ------ Initiating the population ------ ##
function init_pop()
    @inbounds for i in 1:total_population
        apply_age_g(human[i])
    end
end

export init_pop


#Apply Household Size Group to each household

#------------Find House for Each Agent -----------#

@inline function find_a_house()
    avalilable_house = findall(d -> d.AvailableSize>=1 , household)
    return rand(avalilable_house)
end

export find_a_house

#----------put each agent in a house ----------#

function dist_housing()

    #Distributing at least one adult in each housing
    adults = filter(x -> x.AgeGroup >= 4, human)
    for i = 1:total_household
        adults[i].Hid = i
        household[i].AvailableSize -= 1
    end

    left_human = filter(x -> x.Hid == 0, human)

    for j = 1:length(left_human)
        left_id = find_a_house()
            if household[left_id].AvailableSize >=1
            left_human[j].Hid = left_id
            household[left_id].AvailableSize -=1
        elseif household[left_id].AvailableSize == 0
            left_human[j].Hid = left_human[j].Hid
            household[left_id].AvailableSize = household[left_id].AvailableSize
            end
        end
end


function get_gsize(x::Human)
    @match (x.AGC) begin
    1 => 10
    2 => 17
    3 => 14
    4 => 11
    5 => 8
    end
end

function find_repeat(x::Human)
    gsize = get_gsize(x)
    not_same = findall(y->y.Hid!=x.Hid && y.AGC==x.AGC,human)
    x.rep_contact = sample(not_same,gsize,replace = false)
end

function set_repeat()
    for xid in 1:total_population
    x = human[xid]
    find_repeat(x)
    end
end


function set_same_hh()
    for xid in 1:total_population
    x = human[xid]
    x.same_hh = findall(y-> y.Hid == x.Hid && x.id!==y.id, human)
    end
end

function initialize()
    init_human()
    init_household()
    init_pop()
    dist_housing()
    set_same_hh()
    set_repeat()
end
export initialize

function insert_infected(num)
    ## inserts the initial number of active individuls
    l = findall(x -> x.Health == SUSC, human)
    if length(l) > 0
        h = sample(l, num; replace = false)
        @inbounds for i in h
            human[i].Health = ACT
            human[i].swap = ACTT
            human[i].tis = 0
            human[i].exp = 2
            human[i].act_hist = true

        end
    end
end
export insert_infected


reset_params() = reset_params(ModelParameters())

function reset_params(ip::ModelParameters)
    # the p is a global const
    # the ip is an incoming different instance of parameters
    # copy the values from ip to p.
    for x in propertynames(p)
        setfield!(p, x, getfield(ip, x))
    end
end

function get_col_inc(hmatrix, hcol) #here hcol relates to health status
    inth = Int(hcol)
    timevec = zeros(Int64, p.sim_time)
    for r in eachrow(hmatrix)
        idx = findall(x-> r[x-1]!=inth && r[x]==inth, 2:length(r))
        if idx !== nothing
            for i in 1:length(idx)
            timevec[idx[i]] += 1
            end
        end
    end
    return timevec
end

function get_col_prev(hmatrix, hcol)
    inth = Int(hcol)
    timevec = zeros(Int64, p.sim_time)
    for (i,c) in enumerate(eachcol(hmatrix))
        idx = findall(x->x == inth , c)
        if idx !== nothing
            ps = length(c[idx])
            timevec[i] = ps
        end
    end
    return timevec
end

function get_inc_prev(hmatrix) #We get inc and prev of each state

    cols = instances(HEALTH)[1:end - 1] # We don't care about the UNDEF
    inc = zeros(Int64, p.sim_time, length(cols))
    prev = zeros(Int64, p.sim_time, length(cols))
    for i = 1:length(cols)
        inc[:,i] = get_col_inc(hmatrix,cols[i])
        prev[:,i] = get_col_prev(hmatrix,cols[i])
    end
    return inc , prev
end

function collect_df(hmatrix) #Here we are generating a dataframe containing inc and prev of each state at each time

    mdf_inc, mdf_prev = get_inc_prev(hmatrix)
    mdf = hcat(mdf_inc,mdf_prev)
    names_inc = Symbol.(string.((Symbol.(instances(HEALTH)[1:end-1])),"_INC"))
    names_prev = Symbol.(string.((Symbol.(instances(HEALTH)[1:end-1])),"_PREV"))
    names = vcat(names_inc...,names_prev...)
    datf = DataFrame(mdf,names)
    insertcols!(datf,1,:time => 1:p.sim_time)

end

function split_by_age(hmatrix,ags) #generates new hmatrix for each age group
    matx = []
    for i = 1:2
        idx = findall(x -> x ==i , ags)
        push!(matx,view(hmatrix, idx, :))
    end
    return matx
end
export split_by_age



function runsim(simnum, ip::ModelParameters,lat_treat::String, act_treat::String, vaccov::Float64, test::String, trace::Float64, dd::Int64)
    Random.seed!(simnum*777)

    hmatrix = main(ip,lat_treat,act_treat,vaccov,test,trace,dd) #main function that stores results in a matrix

    ags = [x.AGR for x in human]
    all = collect_df(hmatrix)
    spl = split_by_age(hmatrix, ags)
    ag1 = collect_df(spl[1])
    ag2 = collect_df(spl[2])
    insertcols!(all, 1, :sim => simnum); insertcols!(ag1, 1, :sim => simnum); insertcols!(ag2, 1, :sim => simnum);
    return (a=all, g1=ag1, g2=ag2)
end
export runsim

function main(ip::ModelParameters,lat_treat::String, act_treat::String, vaccov::Float64, test::String, trace::Float64, dd::Int64)

reset_params(ip)
get_parameters("Lat_a","Act_a", 0.0, "TST",0.1,4)
hmatrix = zeros(Int16,total_population, p.sim_time)
initialize()
insert_infected(p.initial_infected)
grps = get_grps()

for st = 1:p.sim_time
    kill_infected()
    make_pop_older(st)
    set_scenario(st,lat_treat, act_treat, vaccov, test, trace, dd)
    save_states(st, hmatrix)
    tran_dyn(grps)
    # test_contacts()
    sw = time_update()
end
return hmatrix

end

function set_scenario(sys_time,lat_treat,act_treat,vaccov,test,trace,dd)
    if sys_time == 5773 #Start of year 112 after warmup(100) + calibration(10) + implement at 2020 (2) = 112*52 = 5824
        get_parameters(lat_treat,act_treat,0.0,test,trace,dd)
    end
    if sys_time == 6033 #Start of year (2025) the vaccine becomes available
        get_parameters(lat_treat,act_treat,vaccov,test,trace,dd)
    end
end

function get_parameters(lat_treat::String, act_treat::String, vaccov::Float64, test::String, trace::Float64, dd::Int64)
    
    p.vaccin_cov = vaccov
    p.test_sens = test == "TST" ? 0.70 : 0.90
    p.trace_cov = trace
    p.diagnosis = dd

    if lat_treat == "Lat_a" #9 months of INH - daily or equivalent
            p.α_red = 0.9
            p.lat_treat_time = 39
            p.compliance = 0.86
        elseif lat_treat == "Lat_b"#6 months of INH - daily or equivalent
            p.α_red = 0.01
            p.lat_treat_time = 26
            p.compliance = 0.01
        elseif lat_treat == "Lat_c" #3 months of INH/RMP - daily or equivalent
            p.α_red = 0.64
            p.lat_treat_time = 12
            p.compliance = 0.97
        elseif lat_treat == "Lat_d" #3 months of INH/RPT once weekly or equivalent
            p.α_red = 0.68
            p.lat_treat_time = 12
            p.compliance = 0.96
        elseif lat_treat == "Lat_e" #4 months of RMP daily or equivalent
            p.α_red = 0.63
            p.lat_treat_time = 17
            p.compliance = 0.80
        elseif lat_treat == "Lat_f" #6 months of INH twice weekly or equivalent
            p.α_red = 0.4
            p.lat_treat_time = 26
            p.compliance = 0.55
    end




    if act_treat == "Act_a"
        p.θ = 0.32
        p.act_treat_time = 28
    elseif act_treat == "Act_b" #Higher Success rate 
        p.θ = p.θ
        p.act_treat_time = 16
    end

end


function make_pop_older(sys_time)
    if rem(sys_time,52) == 0
        for xid in 1:total_population
        x = human[xid]
        x.Age += 1
        get_new_AgeGroup(x)
        if x.Age in (5,20,50,65)
            find_repeat(x)
        end

        if rand() < get_death_prob(x)
        x.Health = SUSC
        x.swap = UNDEF
        x.tis = 0
        x.exp = 99999
        x.treated = false
        x.vaccinated = false
        x.act_hist = false
        x.nextweek_meets = 0
        x.will_act = false
        x.HIV = false
        x.Diab = false
        x.renal = false
        apply_age_g(x)
        find_repeat(x)
        end 
        end
    end
    
end

function kill_infected()
    for xid in 1:total_population
        x = human[xid]
        if x.Health == ACT
        if rand()<0.03
        x.Health = SUSC
        x.swap = UNDEF
        x.tis = 0
        x.exp = 99999
        x.treated = false
        x.vaccinated = false
        x.act_hist = false
        x.nextweek_meets = 0
        x.will_act = false
        x.HIV = false
        x.Diab = false
        x.renal = false
        apply_age_g(x)
        find_repeat(x)
            end
        end
    end
end



function tran_dyn(grps)
    totalmet = 0
    totalinf = 0

inf = findall(x-> x.Health == ACT || x.Health== ACTT && x.tis in (1,2) , human )

    beta = p.β

for xid in inf
    x = human[xid]
    
    if x.tis== 2
        beta = (p.β)/2
    end
    
    if x.tis>=3
        x.iso = false
    end


    cnt = x.nextweek_meets

    if length(x.same_hh)!=0 #Splits the number of contacts between hh and community
    cnt_hh = Int(round(cnt*0.85))
    cnt_c = x.iso == true ? 0 : (cnt-cnt_hh)
    else
    cnt_hh = 0
    cnt_c = cnt
    end
    
    cnt_rep = Int(round(cnt_c*0.85,RoundUp)) #This is the contact for repeated contacts
    cnt_left = cnt_c - cnt_rep


    cnt_hh == 0 && @goto c_trans

        if length(x.same_hh)!=0
        meet_hh = rand(x.same_hh,cnt_hh) #sampling contacts from hh
        for j in meet_hh
            y = human[j]
            ycnt = y.nextweek_meets
            ycnt == 0 && continue
            y.nextweek_meets = y.nextweek_meets - 1 #we remove a contact of y
            totalmet += 1

            if y.Health == SUSC && y.swap == UNDEF
                if rand()<beta
                    totalinf += 1
                    y.swap = LAT
                    y.exp = y.tis #The infected individul becomes latent rightaway
                    y.sickfrom = x.id
                    if rand()< (p.test_sens)
                        y.treated = true
                    end
                        
                end
            
            end
        end
    end

    @label c_trans      
    cnt_c == 0 && continue
    cnt_rep == 0 && continue

    @label rep_trans

     meet_rep = rand(x.rep_contact,cnt_rep) #sampling contacts from hh
        for j in meet_rep
            y = human[j]
            ycnt = y.nextweek_meets
            ycnt == 0 && continue
            y.nextweek_meets = y.nextweek_meets - 1 #we remove a contact of y
            totalmet += 1

            if y.Health == SUSC && y.swap == UNDEF
                if rand() <beta
                    totalinf += 1
                    y.swap = LAT
                    y.exp = y.tis #The infected individul becomes latent rightaway
                    y.sickfrom = x.id
                    if rand()<(p.trace_cov)*(p.test_sens)
                        y.treated = true
                    end
                        
                end
            
            end
        end
    


    @label left_tran
    cnt_left == 0 && continue

    gpw = Int.(round.(cm[x.AGC]*cnt_left)) #Based on CM we are devinding the number of contacts
    for (i,g) in enumerate(gpw)
        meet_c = rand(grps[i],g)
        for j in meet_c
            y = human[j]
            ycnt = y.nextweek_meets
            ycnt == 0 && continue
            y.nextweek_meets = y.nextweek_meets - 1 #we remove a contact of y
            totalmet += 1

            if y.Health == SUSC && y.swap == UNDEF
                if rand() <beta
                    totalinf += 1
                    y.swap = LAT
                    y.exp = y.tis #The infected individul becomes latent rightaway
                    y.sickfrom = x.id

                end
            
            end

            
        end
    end
            
end #for end



        return totalmet, totalinf
end #function end
export tran_dyn



# function test_contacts()
#     infector = findall(x-> x.Health == ACTT && x.tis >=3,human)
#     for xid in infector
#         x = human[xid]
#         infectee = findall(y -> y.sickfrom == x.id && y.Hid != x.Hid, human)
#         length(infectee) == 0 && continue
#         tested = rand(infectee,Int(round(p.trace_cov*length(infectee))))
#         for yid in tested
#             if rand()<(p.test_sens)
#                 human[yid].treated = true
#             end
                
#         end
#     end

# end

function get_grps()

    grps = map(x-> findall(y->y.AGC ==x, human), 1:length(Age_braks))
    return grps

end



# function get_latent_dur()
#     latent_p = rand(latent_dist)
#     latent_dur= rand(latent_braks[latent_p])
#     return latent_dur
# end
# export get_latent_dur


function get_latent_dur()
    z = CSV.read("latent.csv")
    latent_p = z[!,1]
    SY = Array{Float64}(undef,130)
    for i = 1:130
        SY[i] = sum(latent_p[(i-1)*4+1:i*4])
    end
    CsumY = cumsum(SY)
    a = rand()
    lat = findfirst(x-> CsumY[x]>=a, 1:130)
    
    if typeof(lat) == Nothing
        lat = 130
    end

    latent_dur = rand((lat-1)*4+1:(lat-1)*4+5)
    return latent_dur
end
export get_latent_dur


function Move_to_Latent(x::Human)

    x.Health = LAT
    x.tis = 0

    lat_dur = get_latent_dur()
    x.swap =  rand() < x.α ? ACT : UNDEF
    x.exp = x.swap == ACT ? lat_dur : 99999

    if x.swap == ACT && x.treated
        x.swap = LATT
        x.exp = rand(1:4) #The person is diagnosed in the first 4 week after infection
    end

end
export Move_to_Latent

function Move_to_Latent_T(x::Human)
    x.Health = LATT
    x.tis = 0
    x.exp = p.lat_treat_time
    x.swap = rand() < p.compliance ? SUSC : DS
    x.α = x.α * (1-p.α_red)
    x.vaccinated = rand()<p.vaccin_cov ? true : false

end
export Move_to_Latent_T

function Move_to_Active(x::Human)
    x.Health = ACT
    x.tis = 0
    x.exp = p.diagnosis #We expect that the person starts treatment after diagnosis time
    x.treated = false
    x.vaccinated = false #we should reset the human
    x.act_hist = true
    x.swap = ACTT

end
export Move_to_Active

function Move_to_Active_T(x::Human)
    x.Health = ACTT
    x.tis = 0
    x.exp = p.act_treat_time
    x.swap = rand()<p.θ ? SUSC : DS #If the treatment is successful the person is SUSC again, o.w, DS
    x.vaccinated = rand()<p.vaccin_cov ? true : false
    x.iso = true

end
export Move_to_Active_T

function Move_to_Dormant(x::Human)
    x.Health = DS
    x.tis = 0
    x.swap = ACT

    if x.act_hist == true && x.vaccinated == true

        x.exp = rand() < 0.46 ? rand(12:84) : 99999 #Reduction in Relapse rate due to vaccination

    elseif x.act_hist == true && x.vaccinated == false

        x.exp = rand(12:84) #No Vaccination, No Reduction in Relapse rate due to vaccination, person gets infected again

    elseif x.act_hist == false && x.vaccinated == true

        x.exp = rand()<x.α *0.55 ? get_latent_dur() : 99999 #Reduction in Activation rate due to vaccination

    elseif x.act_hist == false && x.vaccinated == false

        x.exp = rand()<x.α ? get_latent_dur() : 99999 #No Vaccination, No Reduction in Activation rate

    end
end
export Move_to_Dormant

function Move_to_Susceptible(x::Human)
    #We genererate a totally new individul

    x.Health = SUSC
    x.tis    = 0
    x.swap   = UNDEF 
    x.exp    = 99999
    x.α   = get_alpha(x)
    x.treated = false
    x.vaccinated = false
    x.act_hist = false
    x.nextweek_meets = 0
    x.will_act = false
end
export Move_to_Susceptible

function time_update()
#Counters for different states
lat=0; latt=0; act=0; actt=0; ds=0; susc=0;
for x in human
    x.tis += 1
    if x.tis >= x.exp
        @match Symbol(x.swap) begin
            :LAT => begin Move_to_Latent(x); lat += 1; end
            :LATT => begin Move_to_Latent_T(x); latt += 1; end
            :ACT => begin Move_to_Active(x); act += 1; end
            :ACTT => begin Move_to_Active_T(x); actt += 1; end
            :DS => begin Move_to_Dormant(x); ds +=1 ; end
            :SUSC => begin Move_to_Susceptible(x); susc +=1 end
            _   => error("swap not happened")

        end
    end

    get_nextweek_meets(x)
end
return (lat, latt ,act , actt, ds, susc)
end
export time_update

function save_states(st,hmatrix)
    for i=1:length(human)
        hmatrix[i,st]=Int(human[i].Health)
    end
end
export save_states

function get_nextweek_meets(x::Human)
    ag = x.AGC
    nwm = sum(rand(nbs[ag],7))
    x.nextweek_meets = nwm
    return nwm
end

function contact_matrix()
    # regular contacts, just with 5 age groups.
    #  0-4, 5-19, 20-49, 50-64, 65+
    CM = Array{Array{Float64, 1}, 1}(undef, 5)
    CM[1] = [0.2287, 0.1839, 0.4219, 0.1116, 0.0539]
    CM[2] = [0.0276, 0.5964, 0.2878, 0.0591, 0.0291]
    CM[3] = [0.0376, 0.1454, 0.6253, 0.1423, 0.0494]
    CM[4] = [0.0242, 0.1094, 0.4867, 0.2723, 0.1074]
    CM[5] = [0.0207, 0.1083, 0.4071, 0.2193, 0.2446]
    return CM
end


function negative_binomials()
    ## the means/sd here are calculated using _calc_avgag in Affan's code
    means = [10.21, 16.793, 13.7950, 11.2669, 8.0027]
    sd = [7.65, 11.7201, 10.5045, 9.5935, 6.9638]
    totalbraks = length(means)
    nbinoms = Vector{NegativeBinomial{Float64}}(undef, totalbraks)
    for i = 1:totalbraks
        p = 1 - (sd[i]^2-means[i])/(sd[i]^2)
        r = means[i]^2/(sd[i]^2-means[i])
        nbinoms[i] =  NegativeBinomial(r, p)
    end
    return nbinoms
end
const nbs = negative_binomials()
const cm = contact_matrix()
export negative_binomials, contact_matrix, nbs, cm

function get_death_prob(x::Human)
    @match (x.Age) begin
        0       =>   0.01699 
        1:4     =>   0.00311
        5:9     =>   0.0004
        10:14   =>   0.00453
        15:19   =>   0.01342
        20:24   =>   0.0182
        25:29   =>   0.01455
        30:34   =>   0.01746 
        35:39   =>   0.01011
        40:44   =>   0.01392
        45:49   =>   0.01926
        50:54   =>   0.02582
        55:69   =>   0.04499
        60:64   =>   0.07417
        65:69   =>   0.10719
        70:74   =>   0.1785
        75:79   =>   0.28416
        80:84   =>   0.5065
        85:89   =>   0.51711
        90:100   =>   1
        _ => println("No Age defined")
    end
end

function get_new_AgeGroup(x::Human) #In every year when the person gets older, the agegp changes
    @match (x.Age) begin
        0:4     =>   1
        5:9     =>   2
        10:14   =>   3
        15:19   =>   4
        20:24   =>   5
        25:29   =>   6
        30:34   =>   7
        35:39   =>   8
        40:44   =>   9
        45:49   =>   10
        50:54   =>   11
        55:69   =>   12
        60:64   =>   13
        65:69   =>   14
        70:74   =>   15
        75:79   =>   16
        80:84   =>   17
        85:89   =>   18
        90:99   =>   19

    end
end

#end #module

