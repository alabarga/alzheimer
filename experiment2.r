MethylDataFileName = "./data/meditacion/MeditacionFinalReport_ctrl-bkg.txt"
SampleDescFileName = "./data/meditacion/SampleDescription.csv"

group = "group"
gcase = c("M")   ### Specify the case group index in the sample.txt file (if "concov" is "ON")
gcontrol = c("C")       ### Specify the control group index in the sample.txt file (if "concov" is "ON")

#### Filtramos muestras erroneas o que no peratenecen al estudio
remove_samples = c()

projectName = "meditacion"
#MethylDataFileName = "./data/meditacion/test100.txt"
