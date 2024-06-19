library(dplyr)
url <- "http://smoke.airfire.org/bluesky-daily/output/hysplit-pp/NAM-4km/2014080100/data/fire_locations.csv"
fires <- read.csv(url, stringsAsFactors=FALSE)


# Take the "fires" dataset
#   then filter for type == "WF"
#   then group by state
#   then calculate total area by state
#   then arrange in descending order by total
#   finally, put the result in wildfireAreaByState
fires %>%
  filter(type == "WF") %>%
  group_by(state) %>%
  summarize(total=sum(area, na.rm=TRUE)) %>%
  arrange(desc(total)) ->
  wildfireAreaByState