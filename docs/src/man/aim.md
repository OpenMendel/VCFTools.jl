
# Ancestry Informative Markers

When sample origins are available, one can select SNPs that are most informative at predicting ancestry for your data — the best Ancestry Informative Markers (AIMs).

## Example data

For illustration, we use [chromosome 22 data](http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/) (129 MB) from 1000 genomes project. As show below, this data contains 424147 SNPs and 2504 samples.


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


Next the [population codes](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/) for 1000 genome's samples are explained [here](https://www.internationalgenome.org/category/population/). 

## Preprocess data

To compute AIM markers, we need:
+ Origin (country/continent/population) of every sample stored in a dictionary
+ Complete SNP information. If any sample contains missing genotypes, the given SNP will have p-value of 1. 

The goal is to save the sample origin information into a dictionary where keys are sample IDs and values are origin. Here, the population origin are stored in the file `1000genomes.population.txt` under `VCFTools/test`.


```julia
# cd to test folder
joinpath(pathof(VCFTools), "test")

# read population origin into a dataframe
df = CSV.read("1000genomes.population.txt")

# create dictionary with key = ID, value = population 
sampleID_to_population = Dict{String, String}()
for (id, population) in eachrow(df)
     sampleID_to_population[id] = population
end
sampleID_to_population
```




    Dict{String,String} with 2504 entries:
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

    [32mProgress: 100%|█████████████████████████████████████████| Time: 0:04:58[39m





    424147-element Array{Float64,1}:
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
DataFrame(pvalues=pvals, rank=aim_rank)
```




<table class="data-frame"><thead><tr><th></th><th>pvalues</th><th>rank</th></tr><tr><th></th><th>Float64</th><th>Int64</th></tr></thead><tbody><p>424,147 rows × 2 columns</p><tr><th>1</th><td>1.0</td><td>193453</td></tr><tr><th>2</th><td>1.0</td><td>193454</td></tr><tr><th>3</th><td>1.0</td><td>193455</td></tr><tr><th>4</th><td>1.0</td><td>193456</td></tr><tr><th>5</th><td>1.0</td><td>193457</td></tr><tr><th>6</th><td>1.0</td><td>193458</td></tr><tr><th>7</th><td>1.0</td><td>193459</td></tr><tr><th>8</th><td>1.0</td><td>193460</td></tr><tr><th>9</th><td>1.0</td><td>193461</td></tr><tr><th>10</th><td>1.0</td><td>193462</td></tr><tr><th>11</th><td>5.68468e-120</td><td>39800</td></tr><tr><th>12</th><td>2.20976e-104</td><td>49662</td></tr><tr><th>13</th><td>1.0</td><td>193463</td></tr><tr><th>14</th><td>1.0</td><td>193464</td></tr><tr><th>15</th><td>1.0</td><td>193465</td></tr><tr><th>16</th><td>2.00626e-61</td><td>96498</td></tr><tr><th>17</th><td>1.0</td><td>193466</td></tr><tr><th>18</th><td>1.4409e-78</td><td>73608</td></tr><tr><th>19</th><td>1.0</td><td>193467</td></tr><tr><th>20</th><td>1.0</td><td>193468</td></tr><tr><th>21</th><td>1.17489e-25</td><td>176783</td></tr><tr><th>22</th><td>1.0</td><td>193469</td></tr><tr><th>23</th><td>1.0</td><td>193470</td></tr><tr><th>24</th><td>1.66383e-76</td><td>75835</td></tr><tr><th>25</th><td>1.0</td><td>193471</td></tr><tr><th>26</th><td>1.0</td><td>193472</td></tr><tr><th>27</th><td>7.71708e-107</td><td>47963</td></tr><tr><th>28</th><td>1.18618e-22</td><td>181376</td></tr><tr><th>29</th><td>1.0</td><td>193473</td></tr><tr><th>30</th><td>1.0</td><td>193474</td></tr><tr><th>&vellip;</th><td>&vellip;</td><td>&vellip;</td></tr></tbody></table>


