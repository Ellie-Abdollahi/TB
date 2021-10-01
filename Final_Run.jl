using Distributed
using Base.Filesystem
using DataFrames
using CSV
using Query
using Statistics
using UnicodePlots
using ClusterManagers
using Dates
using DelimitedFiles
## load the packages by TB
using Parameters, Distributions, StatsBase, StaticArrays, Random, Match, DataFrames

addprocs(SlurmManager(500), N=17, topology=:master_worker, exeflags="--project=.")
#addprocs(3)
nsims = 500
# nsims = 3
yr = 150

@everywhere include("TB_main.jl")

myp = ModelParameters()
function final_run(myp::ModelParameters, lat_treat::String, act_treat::String, vaccov::Float64, test::String, trace::Float64, dd::Int64)
    println("starting $nsims simulations")

    # will return 6 dataframes. 1 total, 4 age-specific
    cdr = pmap(1:nsims) do x
            runsim(x, myp,lat_treat,act_treat,vaccov,test, trace, dd) #Here each proc runs 1 sim in cluster by pmap function
    end
    return cdr
end

function yearly_count(col,yr)
cnt = Array{Int64}(undef,yr)
for i in 0:(yr-1)
    cnt[i+1] = round(sum(col[(i*52)+1:52*(i+1)]))
end
return cnt
end


function run_and_save(myp,a,b,c,d,e,f)
cdr = final_run(myp,a,b,c,d,e,f) #This is scenario
allag = vcat([cdr[i].a  for i = 1:nsims]...)
ag1 = vcat([cdr[i].g1 for i = 1:nsims]...)
ag2 = vcat([cdr[i].g2 for i = 1:nsims]...)
#mydfs = Dict("all" => allag)
mydfs = Dict("all" => allag, "ag1" => ag1, "ag2" => ag2)
r1 = Symbol.((:SUSC, :LAT, :LATT, :ACT, :ACTT, :DS), :_INC)
r2 = Symbol.((:SUSC, :LAT, :LATT, :ACT, :ACTT, :DS), :_PREV)
for (k, df) in mydfs
        println("saving dataframe sim level: $k")
        println("saving dataframe yearly level: $k")
        # simulation level, save file per health status, per age group
        for r in vcat(r1..., r2...)
            udf = unstack(df, :time, :sim, r)
            # fn = string("/data/ellie/",string(a),"_",string(b),"_",string(c),"_",string(d),"_",string(e),"/simlevel_", lowercase(string(r)), "_", k, ".csv")
            fn = string("/data/ellie/",string(a),"_",string(b),"_",string(c),"_",string(d),"_",string(e),"_",string(f),"week","/simlevel_", lowercase(string(r)), "_", k, ".csv")
            CSV.write(fn, udf)
            udf = transform(udf[:,2:nsims+1], AsTable(:) => ByRow(mean) => :mean) #Adds mean col to the dataframe
            insertcols!(udf,1,:time=>1:p.sim_time) #Brings the time column back
            ydf = DataFrame()
            insertcols!(ydf,1,:year=>1:yr)
            insertcols!(ydf,2,:mean=>yearly_count(udf.mean,yr))
            for i in 2:nsims+1
                insertcols!(ydf,"$(i-1)"=>yearly_count(udf[:,i],yr))
            end
            # yn = string("/data/ellie/",string(a),"_",string(b),"_",string(c),"_",string(d),"_",string(e),"/yearly_", lowercase(string(r)), "_", k, ".csv")
            yn = string("/data/ellie/",string(a),"_",string(b),"_",string(c),"_",string(d),"_",string(e),"_",string(f),"week","/yearly_", lowercase(string(r)), "_", k, ".csv")
            CSV.write(yn, ydf)
        end
    end
end



run_and_save(myp,"Lat_b","Act_a",0.0,"TST",0.1,4)

# run_and_save(myp,"Lat_a","Act_a",0.0,"TST",0.1,4)
# run_and_save(myp,"Lat_a","Act_a",0.0,"TST",0.3,4)
# run_and_save(myp,"Lat_a","Act_a",0.0,"TST",0.5,4)
# run_and_save(myp,"Lat_a","Act_a",0.0,"TST",0.1,3)
# run_and_save(myp,"Lat_a","Act_a",0.0,"TST",0.3,3)
# run_and_save(myp,"Lat_a","Act_a",0.0,"TST",0.5,3)
# run_and_save(myp,"Lat_a","Act_a",0.0,"TST",0.1,2)
# run_and_save(myp,"Lat_a","Act_a",0.0,"TST",0.3,2)
# run_and_save(myp,"Lat_a","Act_a",0.0,"TST",0.5,2)
# run_and_save(myp,"Lat_a","Act_a",0.0,"TST",0.1,1)
# run_and_save(myp,"Lat_a","Act_a",0.0,"TST",0.3,1)
# run_and_save(myp,"Lat_a","Act_a",0.0,"TST",0.5,1)

# run_and_save(myp,"Lat_a","Act_a",0.0,"IFN",0.1,4)
# run_and_save(myp,"Lat_a","Act_a",0.0,"IFN",0.3,4)
# run_and_save(myp,"Lat_a","Act_a",0.0,"IFN",0.5,4)
# run_and_save(myp,"Lat_a","Act_a",0.0,"IFN",0.1,3)
# run_and_save(myp,"Lat_a","Act_a",0.0,"IFN",0.3,3)
# run_and_save(myp,"Lat_a","Act_a",0.0,"IFN",0.5,3)
# run_and_save(myp,"Lat_a","Act_a",0.0,"IFN",0.1,2)
# run_and_save(myp,"Lat_a","Act_a",0.0,"IFN",0.3,2)
# run_and_save(myp,"Lat_a","Act_a",0.0,"IFN",0.5,2)
# run_and_save(myp,"Lat_a","Act_a",0.0,"IFN",0.1,1)
# run_and_save(myp,"Lat_a","Act_a",0.0,"IFN",0.3,1)
# run_and_save(myp,"Lat_a","Act_a",0.0,"IFN",0.5,1)



# run_and_save(myp,"Lat_b","Act_a",0.0,"TST",0.1,4)
# run_and_save(myp,"Lat_b","Act_a",0.0,"TST",0.3,4)
# run_and_save(myp,"Lat_b","Act_a",0.0,"TST",0.5,4)
# run_and_save(myp,"Lat_b","Act_a",0.0,"TST",0.1,3)
# run_and_save(myp,"Lat_b","Act_a",0.0,"TST",0.3,3)
# run_and_save(myp,"Lat_b","Act_a",0.0,"TST",0.5,3)
# run_and_save(myp,"Lat_b","Act_a",0.0,"TST",0.1,2)
# run_and_save(myp,"Lat_b","Act_a",0.0,"TST",0.3,2)
# run_and_save(myp,"Lat_b","Act_a",0.0,"TST",0.5,2)
# run_and_save(myp,"Lat_b","Act_a",0.0,"TST",0.1,1)
# run_and_save(myp,"Lat_b","Act_a",0.0,"TST",0.3,1)
# run_and_save(myp,"Lat_b","Act_a",0.0,"TST",0.5,1)



# run_and_save(myp,"Lat_c","Act_a",0.0,"TST",0.1,4)
# run_and_save(myp,"Lat_c","Act_a",0.0,"TST",0.3,4)
# run_and_save(myp,"Lat_c","Act_a",0.0,"TST",0.5,4)
# run_and_save(myp,"Lat_c","Act_a",0.0,"TST",0.1,3)
# run_and_save(myp,"Lat_c","Act_a",0.0,"TST",0.3,3)
# run_and_save(myp,"Lat_c","Act_a",0.0,"TST",0.5,3)
# run_and_save(myp,"Lat_c","Act_a",0.0,"TST",0.1,2)
# run_and_save(myp,"Lat_c","Act_a",0.0,"TST",0.3,2)
# run_and_save(myp,"Lat_c","Act_a",0.0,"TST",0.5,2)
# run_and_save(myp,"Lat_c","Act_a",0.0,"TST",0.1,1)
# run_and_save(myp,"Lat_c","Act_a",0.0,"TST",0.3,1)
# run_and_save(myp,"Lat_c","Act_a",0.0,"TST",0.5,1)



# run_and_save(myp,"Lat_d","Act_a",0.0,"TST",0.1,4)
# run_and_save(myp,"Lat_d","Act_a",0.0,"TST",0.3,4)
# run_and_save(myp,"Lat_d","Act_a",0.0,"TST",0.5,4)
# run_and_save(myp,"Lat_d","Act_a",0.0,"TST",0.1,3)
# run_and_save(myp,"Lat_d","Act_a",0.0,"TST",0.3,3)
# run_and_save(myp,"Lat_d","Act_a",0.0,"TST",0.5,3)
# run_and_save(myp,"Lat_d","Act_a",0.0,"TST",0.1,2)
# run_and_save(myp,"Lat_d","Act_a",0.0,"TST",0.3,2)
# run_and_save(myp,"Lat_d","Act_a",0.0,"TST",0.5,2)
# run_and_save(myp,"Lat_d","Act_a",0.0,"TST",0.1,1)
# run_and_save(myp,"Lat_d","Act_a",0.0,"TST",0.3,1)
# run_and_save(myp,"Lat_d","Act_a",0.0,"TST",0.5,1)



# run_and_save(myp,"Lat_e","Act_a",0.0,"TST",0.1,4)
# run_and_save(myp,"Lat_e","Act_a",0.0,"TST",0.3,4)
# run_and_save(myp,"Lat_e","Act_a",0.0,"TST",0.5,4)
# run_and_save(myp,"Lat_e","Act_a",0.0,"TST",0.1,3)
# run_and_save(myp,"Lat_e","Act_a",0.0,"TST",0.3,3)
# run_and_save(myp,"Lat_e","Act_a",0.0,"TST",0.5,3)
# run_and_save(myp,"Lat_e","Act_a",0.0,"TST",0.1,2)
# run_and_save(myp,"Lat_e","Act_a",0.0,"TST",0.3,2)
# run_and_save(myp,"Lat_e","Act_a",0.0,"TST",0.5,2)
# run_and_save(myp,"Lat_e","Act_a",0.0,"TST",0.1,1)
# run_and_save(myp,"Lat_e","Act_a",0.0,"TST",0.3,1)
# run_and_save(myp,"Lat_e","Act_a",0.0,"TST",0.5,1)


# run_and_save(myp,"Lat_a","Act_b",0.0,"TST",0.1,4)
# run_and_save(myp,"Lat_a","Act_b",0.0,"TST",0.3,4)
# run_and_save(myp,"Lat_a","Act_b",0.0,"TST",0.5,4)
# run_and_save(myp,"Lat_a","Act_b",0.0,"TST",0.1,3)
# run_and_save(myp,"Lat_a","Act_b",0.0,"TST",0.3,3)
# run_and_save(myp,"Lat_a","Act_b",0.0,"TST",0.5,3)
# run_and_save(myp,"Lat_a","Act_b",0.0,"TST",0.1,2)
# run_and_save(myp,"Lat_a","Act_b",0.0,"TST",0.3,2)
# run_and_save(myp,"Lat_a","Act_b",0.0,"TST",0.5,2)
# run_and_save(myp,"Lat_a","Act_b",0.0,"TST",0.1,1)
# run_and_save(myp,"Lat_a","Act_b",0.0,"TST",0.3,1)
# run_and_save(myp,"Lat_a","Act_b",0.0,"TST",0.5,1)







