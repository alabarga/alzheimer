MethylDataFileName = "./data/alzheimer/SampleMethFinalReport_total.txt"
SampleDescFileName = "./data/alzheimer/SampleDescription.csv"

group = "group"
gcase = c("EA1", "EA2", "EA3")   ### Specify the case group index in the sample.txt file (if "concov" is "ON")
gcontrol = c("EA_CONTROL")       ### Specify the control group index in the sample.txt file (if "concov" is "ON")

#### Filtramos muestras erroneas o que no pertenecen al estudio
remove_samples = c("46 ER","57 HC")

projectName = 'alzheimer'

