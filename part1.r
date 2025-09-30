setwd("/Users/koo/Desktop")

# 读取莎士比亚文本
a <- scan("shakespeare.txt", what = "character", skip = 83, nlines = 196043 - 83, 
          fileEncoding = "UTF-8")

# The function is split_punct
# Puepose
#  The given punctuation will be separated from the words to make them independent tokens
#  Make sure that no punctuation information is lost and avoid the punctuation marks sticking to the words, which may affect the subsequent processing.
# Inputs
#  words: character vector (sequence of tokens)
#  punct: character vector of punctuation marks to be separated
# Outputs
#  character vector: Token sequences with the same meaning, where punctuation appears as independent elements
# How it works
#  for each punction mark in punct
#   1) Regular expression escaping punctuation
#   2) Find the token that contains this punctuation mark
#   3) For the matched tokens, perform the insertion of "punctuation-free word + punctuation". Tokens that do not match are retained as they are.
#  Update the words step by step until all punctuation is processed
split_punct <- function(words, punct) {
  # Process each punctuation mark
  for (p in punct) {
    # Escape for regular expression matching
    p_esc <- gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", p)
    
    # Find words containing punctuation marks and easy to hand escaped metacharacters
    indices <- grep(p_esc, words, fixed = FALSE)
    
    # Make space for inserting new elements: each hit position will have an additional "punctuation" position.
    new_words <- vector("list", length(words) + length(indices))
    new_idx <- 1
    
    # Traverse all words：If it hits, split it into "pure words + punctuation", otherwise keep it as it is
    for (i in seq_along(words)) {
      if (i %in% indices) {
        # separate pure words and punctuation
        word <- gsub(p_esc, "", words[i])
        # Write pure words and punctuation mark p into the new sequence in sequence
        new_words[[new_idx]] <- word
        new_idx <- new_idx + 1
        new_words[[new_idx]] <- p
        new_idx <- new_idx + 1
      } else {
        new_words[[new_idx]] <- words[i]
        new_idx <- new_idx + 1
      }
    }
    
    # Flatten the list into a character vector and proceed to process the next punctuation mark
    words <- unlist(new_words)
  }
  return(words)
}

# Preprocessed text
# Remove stage directions
# Core idea
#  These explanatory texts are not part of the dialogue content. Incorporating them into the word frequency count would cause deviations.
# Approach
#  Find the indices of all "[" and "]"; for each "[", pair it with the nearest "]" after it and remove the interval
stage_start <- grep("\\[", a)
stage_end <- grep("\\]", a)

stage_indices <- c()
for (i in stage_start) {
  # Find the nearest closing parenthesis
  end <- stage_end[stage_end >= i][1]
  if (!is.na(end)) {
    stage_indices <- c(stage_indices, i:end)
  }
}
# Remove stage directions
if (length(stage_indices) > 0) {
  a <- a[-stage_indices]
}

# Remove all uppercase words and numbers
# Special circumstances: retain "I" and "A"
not_upper <- function(word) {
  if (word %in% c("I", "A")) return(TRUE)  # retain
  word != toupper(word) || grepl("[0-9]", word)  # If it is not all in capital letters, keep it as it is or if it contains numbers, keep it as it is.
}
keep <- sapply(a, not_upper)
a <- a[keep]

# Remove the hyphen
a <- gsub("-", "", a, fixed = TRUE)
a <- gsub("—", "", a, fixed = TRUE)

# Call the function(split_punct) to separate punctuation
punct_marks <- c(",", ".", ";", ":", "!", "?")
a <- split_punct(a, punct_marks)

# Convert to lowercase (except for "i" and "a")
to_lower_special <- function(word) {
  if (word %in% c("I", "A")) return(word)
  tolower(word)
}
a <- sapply(a, to_lower_special)
names(a) <- NULL  # Remove the name attribute

# Create common word vectors
# Build a unique vocabulary list
b <- unique(a)

# Map each word in a to its position in b
a_indices <- match(a, b)

# Count word frequencies based on the index statistics
word_counts <- tabulate(a_indices)

# Select the 1,000 most common high-frequency words
k <- 1000
r <- rank(-word_counts, ties.method = "first")   # The larger the negative sign, the higher the rank
top_k_indices <- which(r <= k)                   # The top k names
b_common <- b[top_k_indices]

