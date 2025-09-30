#Yunhan Zhang s2176155. Xiyu Wu s2799746. Tianyu Wang s2794991
#Yunhan did question 6 and 7. Xiyu did question 4 and 5. Tianyu did question 8 and 9.
#Each member did roughly same amount work.

setwd("C:\\Users\\24136\\Desktop\\UoE-postgra\\ESP\\ESP_37")
#setwd("/Users/koo/Desktop")

# Read the file
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


# Convert the cleaned text 'a' into a vector of integer tokens.
# Each word in 'a' is matched to its index in b_common (the top ~1000 most frequent words).
# If a word is not in b_common, match() returns NA for that position.
tokens <- match(a, b_common)

# Set the maximum order (mlag). 
# For mlag = 4 means using the previous 4 words to predict the 5th.
mlag <- 4  

# Total number of tokens in the text
n <- length(tokens)

# Pre-allocate an (n - mlag) x (mlag + 1) matrix M.
M <- matrix(NA, nrow = n - mlag, ncol = mlag + 1)

# Fill M column by column using sliding windows of the token vector.
# - Column 1 = tokens[1:(n-mlag)]
# - Column 2 = tokens[2:(n-mlag+1)]
# - ...
# - Column (mlag+1) = tokens[(mlag+1):n]
for (i in 1:(mlag + 1)) {
  M[, i] <- tokens[i:(n - mlag + i - 1)]
}

# Removing rows containing NA ensures that every sequence in M consists only of "common words",
M <- M[complete.cases(M), ]

# ================================================================
# Next-word prediction function
# ------------------------------------------------
# Purpose:
# Given the recent tokens, this function predicts
# the next token by searching for matching contexts in M.
#
# Inputs:
# - key: integer vector of recent tokens (the current context)
# - M: the sequence matrix built above (contexts + next words)
# - M1: the full vector of tokens for the entire text
# - w: optional vector of mixture weights for each order (default equal)
#
# Output:
# - An integer representing the predicted next token, chosen at random
#   according to the estimated probability distribution.
#
# How it works :
# - Try matching the last i tokens of 'key' against M for i = 1..mlag.
# - Collect all possible next tokens that followed such contexts in Shakespeare.
# - Build a probability distribution across those tokens, weighted by w[i].
# - Sample one token from this distribution.
# - If no matches are found, fall back to unigram frequencies from M1.
# ================================================================
next_word <- function(key, M, M1, w = rep(1, ncol(M) - 1)) {
  u <- integer(0)  # candidate next tokens
  p <- numeric(0)  # probabilities
  
  # Check all possible context lengths (orders), from 1 up to mlag.
  max_i <- min(length(key), ncol(M) - 1)
  for (i in 1:max_i) {
    # Extract the last i tokens from the current key
    current_key <- key[(length(key) - i + 1):length(key)]
    
    # Match against the last i context columns in M
    mc <- (ncol(M) - 1) - i + 1
    # Find the rows of M that match key
    ii <- colSums(!(t(M[, mc:(ncol(M) - 1), drop = FALSE]) == current_key))
    row_match <- (ii == 0)   # rows where context matches exactly
    
    if (any(row_match)) {
      # Collect the next tokens that actually followed this context in Shakespeare
      next_words <- M[row_match, ncol(M)]
      
      # Share the mixture weight for this order equally among all matches
      prob_each <- w[i] / length(next_words)
      
      # Append candidates and their probabilities
      u <- c(u, next_words)
      p <- c(p, rep(prob_each, length(next_words)))
    }
  }
  
  # If no matches were found at any order, fall back to unigram distribution
  if (length(u) == 0) {
    K <- max(M1, na.rm = TRUE)
    freq <- tabulate(M1[!is.na(M1)], nbins = K)   # word counts
    prob <- if (sum(freq) > 0) freq / sum(freq) else rep(1 / K, K)
    return(sample.int(K, size = 1, prob = prob))
  }
  
  # Normalize probabilities and sample one next token
  p <- p / sum(p)
  sample(u, size = 1, prob = p)
}