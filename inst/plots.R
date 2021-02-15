library(kwdb)
library(kwdemog)
library(dplyr)
library(ggplot2)
library(ggsidekick)
library(viridis)

data(kwdata)

# expand data
expanded_data = expand(kwdata)

# create age class
expanded_data$class = NA
expanded_data$class[which(expanded_data$age < 10)] = "Juvenile"
expanded_data$class[which(expanded_data$age %in% seq(10,40) & expanded_data$sexF1M2 == 1)] = "Repro.F"
expanded_data$class[which(expanded_data$age > 40 & expanded_data$sexF1M2 == 1)] = "Old.F"
expanded_data$class[which(expanded_data$age %in% seq(10,22) & expanded_data$sexF1M2 == 2)] = "Young.M"
expanded_data$class[which(expanded_data$age > 22 & expanded_data$sexF1M2 == 2)] = "Old.M"

expanded_data$class = factor(expanded_data$class)
expanded_data$class = factor(expanded_data$class,levels(expanded_data$class)[c(1,4,2,5,3)])

srkw_deaths = dplyr::filter(expanded_data, population == "SRKW") %>%
  dplyr::filter(alive == 0)


srkw_deaths$period = NA
srkw_deaths$period[which(srkw_deaths$year %in% 1995:2001)] = "1995-01"
srkw_deaths$period[which(srkw_deaths$year %in% 2011:2018)] = "2011-18"

p1 = filter(srkw_deaths, is.na(period)==FALSE) %>%
  ggplot(aes(class), group=period) + geom_bar(aes(fill=period), position='dodge') + theme_sleek() +
  scale_fill_viridis(discrete = TRUE,begin=0.1, end=0.5) +
  scale_color_viridis(discrete = TRUE, begin=0.1, end=0.5) +
  xlab("Age-sex class") + ylab("Count") +
  ggtitle("SRKW Mortalities (pre-ESA listing, and since 2011)") +
  theme(text = element_text(size=20))
pdf("Mortalities_by_era.pdf")
print(p1)
dev.off()


srkw_alive = dplyr::filter(expanded_data, population == "SRKW") %>% filter(alive==1) %>%
  group_by(year, class) %>%
  summarize(n = n()) %>%
  ggplot(aes(year, n, group=class, color = class)) + geom_line(size=2) +
  scale_fill_viridis(discrete = TRUE,begin=0, end=0.7) +
  scale_color_viridis(discrete = TRUE, begin=0, end=0.7) +
  ylab("Count") + xlab("Year") +
  ggtitle("SRKW age and sex structure") +
  theme_sleek() + theme(text = element_text(size=20))
pdf("SRKW_age_structure.pdf", width=10,height= 7)
print(srkw_alive)
dev.off()

nrkw_alive = dplyr::filter(expanded_data, population == "NRKW") %>% filter(alive==1) %>%
  group_by(year, class) %>%
  summarize(n = n()) %>%
  ggplot(aes(year, n, group=class, color = class)) + geom_line(size=2) +
  scale_fill_viridis(discrete = TRUE,begin=0, end=0.7) +
  scale_color_viridis(discrete = TRUE, begin=0, end=0.7) +
  ylab("Count") +
  ggtitle("NRKW age and sex structure") +
  theme(text = element_text(size=20)) + theme_sleek()


srkw_births = dplyr::filter(kwdata) %>% filter(sexF1M2 !=0) %>%
  filter(birth > 1976) %>%
  mutate(sex = sexF1M2 - 1) %>%
  ggplot(aes(birth, sex)) + geom_point(col="dark blue", alpha=0.5, size=4) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  ylab("Sex at birth (0 = female)") +theme_sleek() +
  theme(text = element_text(size=20)) + xlab("") +
  geom_hline(aes(yintercept=0.5), col="red") + facet_wrap(~ population)
pdf("SRKW_sex_ratio.pdf", height=7, width=10)
print(srkw_births)
dev.off()


srkw_alive = dplyr::filter(expanded_data, population == "SRKW") %>% filter(alive==1) %>% filter(sexF1M2 == 1) %>%
  filter(year ==2017) %>% filter(age %in% seq(10,42))
srkw_alive$last_birth = NA
for(i in 1:nrow(srkw_alive)) {
  this_animal = filter(expanded_data, animal == srkw_alive$animal[i])
  births = which(this_animal$gave_birth==1)
  if(length(births) > 0) {
    srkw_alive$last_birth[i] = max(this_animal$age[births])
  }
}
srkw_alive = srkw_alive[,c("animal","pod","age","last_birth")]
srkw_alive$diff = srkw_alive$age - srkw_alive$last_births
srkw_alive = arrange(srkw_alive, pod, -age)
srkw_alive$pod[which(srkw_alive$pod == "J001")] = "J"
srkw_alive$pod[which(srkw_alive$pod == "K001")] = "K"
srkw_alive$pod[which(srkw_alive$pod == "L001")] = "L"
srkw_alive = select(srkw_alive, -last_birth)



