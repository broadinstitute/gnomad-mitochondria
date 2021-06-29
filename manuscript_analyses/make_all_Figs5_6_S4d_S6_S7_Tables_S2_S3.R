#do analyses and make figures related to patterns of variation, and pathogenic variants
dir.create("figures")
dir.create("tables")

#main text figures
source("collate_Figure5.R", echo = TRUE) 
source("make_Fig6.R", echo = TRUE) 

#supplementary figures
source("make_FigS4d.R", echo = TRUE) 
source("make_FigS6.R", echo = TRUE) 
source("make_FigS7.R", echo = TRUE)
#options(bitmapType = 'cairo', device = 'png')
ggsave("figures/FigureS7.png", width = 15, height = 13) #due to a viewport error it isn't saving the figure within the script, so providing command here

#supplementary tables
source("make_TableS2.R", echo = TRUE) 
source("make_TableS3.R", echo = TRUE) 





