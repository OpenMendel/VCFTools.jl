
# Ancestry Informative Markers

When sample origins are available, one can select SNPs that are most informative at predicting ancestry for your data — the best Ancestry Informative Markers (AIMs).

## Example data

For illustration, we use [chromosome 22 data](http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/) (129 MB) from 1000 genomes project. As shown below, this data contains 424147 SNPs and 2504 samples.


```julia
using Revise
using VCFTools
using Random
using CSV
using DataFrames
using StatsBase
```


```julia
# download data
vcffile = "chr22.1kg.phase3.v5a.vcf.gz"
isfile(vcffile) || 
    download("http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/chr22.1kg.phase3.v5a.vcf.gz", 
    joinpath(pwd(), vcffile))

# compute simple summary statistics
@show nrecords(vcffile)
@show nsamples(vcffile);
```

    nrecords(vcffile) = 424147
    nsamples(vcffile) = 2504


Note the [population codes](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/) for 1000 genome's samples are explained [here](https://www.internationalgenome.org/category/population/). 

## Preprocess data

To compute AIM markers, we need:
+ Origin (country/continent/population) of every sample stored in a dictionary
+ Complete SNP information. If any sample contains missing genotypes, the given SNP will have p-value of 1. 

The goal is to save the sample origin information into a dictionary where keys are sample IDs and values are origin. Here, the population origin are stored in the file `1000genomes.population.txt` under `VCFTools/test`.


```julia
# cd to test folder
cd(normpath(VCFTools.datadir()))

# read population origin into a dataframe
df = CSV.read("1000genomes.population.txt", DataFrame)

# create dictionary with key = ID, value = population 
sampleID_to_population = Dict{String, String}()
for (id, population) in eachrow(df)
     sampleID_to_population[id] = population
end
sampleID_to_population
```




    Dict{String, String} with 2504 entries:
      "HG01791" => "GBR"
      "HG02736" => "PJL"
      "HG00182" => "FIN"
      "HG03914" => "BEB"
      "HG00149" => "GBR"
      "NA12156" => "CEU"
      "HG02642" => "GWD"
      "HG02851" => "GWD"
      "NA19835" => "ASW"
      "NA19019" => "LWK"
      "HG01131" => "CLM"
      "HG03578" => "MSL"
      "NA18550" => "CHB"
      "HG02401" => "CDX"
      "HG01350" => "CLM"
      "HG03973" => "ITU"
      "NA07000" => "CEU"
      "HG01709" => "IBS"
      "HG01395" => "PUR"
      "HG01980" => "PEL"
      "HG01979" => "PEL"
      "HG01122" => "CLM"
      "HG03869" => "ITU"
      "HG03729" => "ITU"
      "NA19920" => "ASW"
      ⋮         => ⋮



## Run AIM selection

A p-value will be computed for each SNP. Smaller p-values indicate more ancestry informative.


```julia
pvals = VCFTools.aim_select(vcffile, sampleID_to_population)
```

    [32mProgress: 100%|█████████████████████████████████████████| Time: 0:05:45[39m





    424147-element Vector{Float64}:
     1.0
     1.0
     1.0
     1.0
     1.0
     1.0
     1.0
     1.0
     1.0
     1.0
     5.684678238417218e-120
     2.2097628180404582e-104
     1.0
     ⋮
     1.0
     1.0
     3.918650824065397e-73
     2.2771505483860234e-77
     1.0
     1.8351937329381248e-92
     1.0
     1.0
     1.0
     3.0769656913253715e-14
     1.0
     1.0



We can rank these p-values


```julia
aim_rank = ordinalrank(pvals)
df = DataFrame(pvalues=pvals, rank=aim_rank)
@show df[1:20, :]; # list 20 rows
```

    df[1:20, :] = 20×2 DataFrame
     Row │ pvalues       rank
         │ Float64       Int64
    ─────┼──────────────────────
       1 │ 1.0           193453
       2 │ 1.0           193454
       3 │ 1.0           193455
       4 │ 1.0           193456
       5 │ 1.0           193457
       6 │ 1.0           193458
       7 │ 1.0           193459
       8 │ 1.0           193460
       9 │ 1.0           193461
      10 │ 1.0           193462
      11 │ 5.68468e-120   39800
      12 │ 2.20976e-104   49662
      13 │ 1.0           193463
      14 │ 1.0           193464
      15 │ 1.0           193465
      16 │ 2.00626e-61    96498
      17 │ 1.0           193466
      18 │ 1.4409e-78     73608
      19 │ 1.0           193467
      20 │ 1.0           193468

