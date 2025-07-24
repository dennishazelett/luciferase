using Turing, Distributions, StatsPlots, CSV, DataFrames, CategoricalArrays


lucdat = CSV.read("luciferase_assays.csv",DataFrame)

lucdat.plate = categorical(lucdat.plate)
lucdat.enhancer = categorical(lucdat.enhancer)
lucdat.plasmid = categorical(lucdat.plasmid)



@model function lucmodel(firefly, renilla, col, row, plate, trans, enh, plas, allele, batch, DHT)
    foo ~ Normal(0.0,1.0)

end

ourmodel = lucmodel(lucdat.firefly,lucdat.renilla,lucdat.column, lucdat.row,
    lucdat.plate, lucdat.transfection, lucdat.enhancer, lucdat.plasmid, 
    lucdat.allele, lucdat.batch, lucdat.DHT)

sam = Turing.sample(ourmodel,NUTS(300,0.85,max_depth=10),MCMCThreads(),500,3)

