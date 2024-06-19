# biopsia liquida
# comparacion controles alzheimer vs meditacion

projectName = 'liquid_sanos'

MethylDataFileName = "./data/liquid/BiopsiaLiquidaMethylationReport_noheader.csv"
SampleDescFileName = "./data/liquid/SampleDescription.csv"

group = "Source"
gcase = c("Brain")   ### Specify the case group index in the sample.txt file (if "concov" is "ON")
gcontrol = c("Blood")       ### Specify the control group index in the sample.txt file (if "concov" is "ON")

#### Filtramos muestras erroneas o que no pertenecen al estudio
remove_samples = c("46 ER","57 HC", as.character(sampleDescriptions$SampleID[(sampleDescriptions$Status == "Target") & (sampleDescriptions$Source=="Brain")]))

