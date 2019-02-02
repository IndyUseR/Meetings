
# Load libraries
library(tidyverse)
library(gt)

# load survey data
surveydf = read_csv("2019_01_January/IndyUseRSurvey.csv")

# gather all questions/answers in tall format
surveydf = surveydf %>% select(-Timestamp) %>% mutate(id = seq(1, nrow(surveydf)))
# split responses where more than one response is given to a question
surveydf_g = surveydf %>% gather(question, answer, -id)
surveydf_g = surveydf_g %>% mutate(anssplit = strsplit(answer,";")) %>% unnest(anssplit)

# question list
questions = surveydf_g %>% distinct(question)
questions = questions$question

getbarplot = function(quest) {
  p = ggplot(surveydf_g %>% filter(question == quest) %>% count(anssplit)) + 
    geom_col(aes(x = fct_reorder(anssplit, n), y = n), width = 0.5, fill = "lightblue") + 
    xlab("") + ylab("# responses") + 
    coord_flip() + theme_bw() + ggtitle(quest)
  return(p)
}

barplotq = c(1, 2, 3, 4, 6, 7, 8, 10, 12)
for(i in barplotq) {
  p = getbarplot(questions[i])
  print(p)
}

getOpenCnt = function(quest) {
   resp = surveydf_g %>% filter(question == quest) %>% 
    filter(!is.na(anssplit)) %>% count(anssplit) %>% pull(anssplit) 
   return(resp)
}

openrespq = surveydf_g %>% filter(question %in% questions[c(5, 9, 11, 13, 14)]) %>% 
              filter(!is.na(anssplit)) %>% select(question, anssplit) %>% 
              rename(response = anssplit)
openrespq %>% mutate(dummy = "         ") %>% gt(rowname_col = "dummy", groupname_col = "question")


