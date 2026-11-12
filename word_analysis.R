library(tidyverse)
library(readxl)
library(tidytext)

raw <- read_xlsx("words_for_analysis.xlsx") |>
  mutate(full_text = paste(full_text, full_text_contd)) |>
  select(!full_text_contd) |>
  mutate(id = paste(author, year))

tb_stop_words <- tibble(word = c("mg", "kg", "dose", "grade"))
words_of_interest <- tibble(
  word = c("safe", "safety", "tolerable", "tolerability", "tolerated",
           "acceptable", "manageable", "feasible")
)

# Abstract
abs1 <- raw |>
  select(!c(full_text, author, year, discussion)) |>
  unnest_tokens(word, abstract) |>
  group_by(id) |>
  add_count(word, sort = T) |>
  distinct(id, word, .keep_all = T) |>
  inner_join(words_of_interest, by = join_by(word))

# Just looking at the bigram well tolerated
abs2 <- raw |>
  select(!c(full_text, author, year, discussion)) |>
  unnest_tokens(word, abstract, token = "ngrams", n = 2) |>
  group_by(id) |>
  filter(word %in% c("well tolerated", 
                     "favourable toxicity",
                     "favorable toxicity",
                     "favourable safety",
                     "favorable safety")) |>
  add_count(word, sort = T) |>
  distinct(id, word, .keep_all = T)

abstract_for_manuscript <- 
  bind_rows(abs1, abs2) |>
  rename(Regimen = Studies, Study = id, 
         Word = word, N = n) |> 
  filter(!is.na(N)) |> 
  mutate(Word = str_to_sentence(Word)) |> 
  ungroup() |> 
  group_by(Regimen, Study) |> 
  arrange(desc(N), .by_group = T) |>
  mutate(year = str_extract(Study, "\\d.*")) |> 
  mutate(Study = str_c(
    str_extract(Study, "\\w.*[^[:digit:]|\\s]"), " et al., ", year)) |> 
  arrange(desc(Regimen), year, Study, Word) |> 
  select(!year)

abstract_for_manuscript |> 
  print(n = Inf)

abstract_for_manuscript |> 
  ungroup() |> 
  count(Word)

# Discussion
discussion_words <- raw |>
  select(!c(abstract, full_text, author, year)) |>
  unnest_tokens(word, discussion) |>
  anti_join(stop_words, by = join_by(word)) |>
  anti_join(tb_stop_words, by = join_by(word)) |>
  filter(!str_detect(word, "\\d")) |>
  group_by(id) |>
  add_count(word, sort = T) |>
  distinct(id, word, .keep_all = T) |>
  # bind_tf_idf(word, id, n) |>
  # arrange(desc(tf_idf)) |>
  right_join(words_of_interest, by = join_by(word)) |>
  arrange(id) |>
  print(n = Inf)

discussion2 <- raw |>
  select(!c(author, year, abstract, full_text)) |>
  unnest_tokens(word, discussion, token = "ngrams", n = 2) |>
  group_by(id) |>
  filter(word %in% c("well tolerated", 
                     "favourable toxicity",
                     "favorable toxicity",
                     "favourable safety",
                     "favorable safety")) |>
  add_count(word, sort = T) |>
  distinct(id, word, .keep_all = T)

discussion_for_manuscript <- 
  bind_rows(discussion_words, discussion2) |>
  rename(Regimen = Studies, Study = id, 
         Word = word, N = n) |> 
  filter(!is.na(N)) |> 
  mutate(Word = str_to_sentence(Word)) |> 
  ungroup() |> 
  group_by(Regimen, Study) |> 
  arrange(desc(N), .by_group = T) |> 
  mutate(year = str_extract(Study, "\\d.*")) |> 
  mutate(Study = str_c(
    str_extract(Study, "\\w.*[^[:digit:]|\\s]"), " et al., ", year)) |> 
  arrange(desc(Regimen), year, Study, Word) |> 
  select(!year)

discussion_for_manuscript |> 
  print(n = Inf)

discussion_for_manuscript |> 
  ungroup() |> 
  count(Word)

write_csv(discussion_for_manuscript, "Discussion Words.csv")
write_csv(abstract_for_manuscript, "Abstract Words.csv")
